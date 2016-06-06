#include <config.h>

#include "msfem_solver.hh"

#include <sstream>
#include <assert.h>
#include <boost/assert.hpp>

#include <dune/common/exceptions.hh>
#include <dune/common/timer.hh>

#include <dune/multiscale/msfem/coarse_scale_operator.hh>
#include <dune/multiscale/msfem/localproblems/localgridsearch.hh>
#include <dune/multiscale/tools/misc.hh>
#include <dune/multiscale/common/traits.hh>
#include <dune/multiscale/common/heterogenous.hh>
#include <dune/multiscale/problems/selector.hh>
#include <dune/multiscale/msfem/msfem_traits.hh>
#include <dune/multiscale/msfem/localproblems/localgridlist.hh>
#include <dune/multiscale/tools/misc/outputparameter.hh>
#include <dune/multiscale/common/df_io.hh>
#include <dune/multiscale/msfem/localproblems/localsolutionmanager.hh>
#include <dune/multiscale/msfem/localproblems/localproblemsolver.hh>

#include <dune/xt/common/logging.hh>
#include <dune/xt/common/configuration.hh>
#include <dune/xt/common/timings.hh>
#include <dune/xt/common/ranges.hh>

#include <dune/gdt/operators/prolongations.hh>
#include <dune/gdt/spaces/cg.hh>
#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/functionals/l2.hh>
#include <dune/gdt/spaces/constraints.hh>
#include <dune/gdt/assembler/system.hh>
#include <dune/gdt/operators/projections.hh>

#include <dune/multiscale/msfem/localsolution_proxy.hh>

namespace Dune {
namespace Multiscale {

void Elliptic_MsFEM_Solver::identify_fine_scale_part(const Problem::ProblemContainer& problem,
                                                     LocalGridList& localgrid_list,
                                                     const CommonTraits::DiscreteFunctionType& coarse_msfem_solution,
                                                     const CommonTraits::SpaceType& coarse_space,
                                                     std::unique_ptr<LocalsolutionProxy>& msfem_solution) const {
  Dune::XT::Common::ScopedTiming st("msfem.idFine");

  auto& coarse_indexset = coarse_space.grid_view().grid().leafIndexSet();
  const bool is_simplex_grid = DSG::is_simplex_grid(coarse_space);

  LocalsolutionProxy::CorrectionsMapType local_corrections;
  const auto interior = coarse_space.grid_view().grid().leafGridView<InteriorBorder_Partition>();
  for (const auto& coarse_entity : Dune::elements(interior)) {
    LocalSolutionManager localSolutionManager(coarse_space, coarse_entity, localgrid_list);
    localSolutionManager.load();
    auto& localSolutions = localSolutionManager.getLocalSolutions();
    const auto coarse_index = coarse_indexset.index(coarse_entity);
    local_corrections[coarse_index] = Dune::XT::Common::make_unique<MsFEMTraits::LocalGridDiscreteFunctionType>(
        localSolutionManager.space(), "correction");

    auto& local_correction = *local_corrections[coarse_index];
    local_correction.vector() *= 0;
    const auto coarseSolutionLF = coarse_msfem_solution.local_discrete_function(coarse_entity);

    if (is_simplex_grid) {
      BOOST_ASSERT_MSG(localSolutions.size() == CommonTraits::world_dim,
                       "We should have dim local solutions per coarse element on triangular meshes!");

      MsFEMTraits::LocalGridDiscreteFunctionType::JacobianRangeType grad_coarse_msfem_on_entity;
      // We only need the gradient of the coarse scale part on the element, which is a constant.
      coarseSolutionLF->jacobian(coarse_entity.geometry().center(), grad_coarse_msfem_on_entity);

      // get the coarse gradient on T, multiply it with the local correctors and sum it up.
      for (int spaceDimension = 0; spaceDimension < CommonTraits::world_dim; ++spaceDimension) {
        localSolutions[spaceDimension]->vector() *= grad_coarse_msfem_on_entity[0][spaceDimension];
        local_correction.vector() += localSolutions[spaceDimension]->vector();
      }
      BOOST_ASSERT_MSG(false, "no adding of the boundary correctors??");
    } else {
      //! @warning At this point, we assume to have the same types of elements in the coarse and fine grid!
      //            BOOST_ASSERT_MSG(
      //                static_cast<long long>(localSolutions.size() - localSolManager.numBoundaryCorrectors()) ==
      //                    static_cast<long long>(coarseSolutionLF.size()),
      //                "The current implementation relies on having thesame types of elements on coarse and fine
      // level!");
      for (std::size_t dof = 0; dof < coarseSolutionLF->vector().size(); ++dof) {
        localSolutions[dof]->vector() *= coarseSolutionLF->vector().get(dof);
        local_correction.vector() += localSolutions[dof]->vector();
      }

      // oversampling : restrict the local correctors to the element T
      // ie set all dofs not "covered" by the coarse cell to 0
      // also adds lg-prolongation of coarse_solution to local_correction
      const auto cut_overlay = problem.config().get("msfem.oversampling_layers", 0) > 0;
      for (const auto& local_entity : Dune::elements(localSolutionManager.space().grid_view())) {
        const auto& lg_points = localSolutionManager.space().lagrange_points(local_entity);
        const auto& reference_element = DSG::reference_element(coarse_entity);
        const auto& coarse_geometry = coarse_entity.geometry();
        auto entity_local_correction = local_correction.local_discrete_function(local_entity);

        for (const auto lg_i : Dune::XT::Common::value_range(int(lg_points.size()))) {
          const auto local_point = lg_points[lg_i];
          const auto global_lg_point = local_entity.geometry().global(local_point);
          const auto local_coarse_point = coarse_geometry.local(global_lg_point);
          auto& vec = entity_local_correction->vector();
          if (cut_overlay) {
            const bool covered = reference_element.checkInside(local_coarse_point);
            vec.set(lg_i, covered ? vec.get(lg_i) : 0);
          }
        }
      }

      // add dirichlet corrector
      //      local_correction.vector() += localSolutions[coarseSolutionLF->vector().size() + 1]->vector();
      // substract neumann corrector
      local_correction.vector() -= localSolutions[coarseSolutionLF->vector().size()]->vector();

      if (problem.config().get("msfem.local_corrections_vtk_output", false)) {
        const std::string name = (boost::format("local_correction_%d_") % coarse_index).str();
        Dune::Multiscale::OutputParameters outputparam(problem.config().get("global.datadir", "data"));
        outputparam.set_prefix(name);
        local_correction.visualize(outputparam.fullpath(local_correction.name()));
      }
    }
    localSolutions.clear();
  }

  MS_LOG_INFO_0 << "Dirichlet correctors are broken and disabled\n";
  msfem_solution =
      Dune::XT::Common::make_unique<LocalsolutionProxy>(std::move(local_corrections), coarse_space, localgrid_list);
}

void Elliptic_MsFEM_Solver::apply(DMP::ProblemContainer& problem, const CommonTraits::SpaceType& coarse_space,
                                  std::unique_ptr<LocalsolutionProxy>& solution, LocalGridList& localgrid_list) const {
  Dune::XT::Common::ScopedTiming st("msfem.Elliptic_MsFEM_Solver.apply");
  const auto clearGuard = DiscreteFunctionIO::clear_guard();

  CommonTraits::DiscreteFunctionType coarse_msfem_solution(coarse_space, "Coarse Part MsFEM Solution");
  coarse_msfem_solution.vector() *= 0;

  //! Solutions are kept in-memory via DiscreteFunctionIO::MemoryBackend by LocalsolutionManagers
  LocalProblemSolver(problem, coarse_space, localgrid_list).solve_for_all_cells();

  CoarseScaleOperator elliptic_msfem_op(problem, coarse_space, localgrid_list);
  elliptic_msfem_op.apply_inverse(coarse_msfem_solution);

  //! identify fine scale part of MsFEM solution (including the projection!)
  identify_fine_scale_part(problem, localgrid_list, coarse_msfem_solution, coarse_space, solution);
  solution->add(coarse_msfem_solution);
}

} // namespace Multiscale {
} // namespace Dune {
