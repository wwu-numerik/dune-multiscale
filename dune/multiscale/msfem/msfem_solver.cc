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
#include <dune/multiscale/problems/selector.hh>
#include <dune/multiscale/msfem/msfem_traits.hh>
#include <dune/multiscale/msfem/localproblems/localgridlist.hh>
#include <dune/multiscale/tools/misc/outputparameter.hh>
#include <dune/multiscale/tools/discretefunctionwriter.hh>
#include <dune/multiscale/msfem/localproblems/localsolutionmanager.hh>
#include <dune/multiscale/msfem/localproblems/localproblemsolver.hh>

#include <dune/stuff/fem/functions/checks.hh>
#include <dune/stuff/common/logging.hh>
#include <dune/stuff/common/configuration.hh>
#include <dune/stuff/common/profiler.hh>
#include <dune/stuff/common/ranges.hh>
#include <dune/stuff/discretefunction/projection/heterogenous.hh>

#include <dune/gdt/operators/prolongations.hh>
#include <dune/gdt/spaces/continuouslagrange.hh>
#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/functionals/l2.hh>
#include <dune/gdt/spaces/constraints.hh>
#include <dune/gdt/assembler/system.hh>
#include <dune/gdt/operators/projections.hh>

#include <dune/multiscale/msfem/localsolution_proxy.hh>

namespace Dune {
namespace Multiscale {

void Elliptic_MsFEM_Solver::identify_fine_scale_part(LocalGridList& subgrid_list,
                                                     const DiscreteFunctionType& coarse_msfem_solution,
                                                     std::unique_ptr<LocalsolutionProxy>& msfem_solution) const {
  DSC::Profiler::ScopedTiming st("msfem.idFine");

  const auto& coarse_space = coarse_msfem_solution.space();
  auto& coarse_indexset = coarse_space.grid_view()->grid().leafIndexSet();
  const bool is_simplex_grid = DSG::is_simplex_grid(coarse_space);

  LocalsolutionProxy::CorrectionsMapType local_corrections;
  for (auto& coarse_entity : DSC::viewRange(*coarse_space.grid_view())) {
    LocalSolutionManager localSolManager(coarse_space, coarse_entity, subgrid_list);
    localSolManager.load();
    auto& localSolutions = localSolManager.getLocalSolutions();
    const auto coarse_index = coarse_indexset.index(coarse_entity);
    local_corrections[coarse_index] =
        DSC::make_unique<MsFEMTraits::LocalGridDiscreteFunctionType>(localSolManager.space(), " ");

    auto& local_correction = *local_corrections[coarse_index];
    local_correction.vector() *= 0;
    const auto coarseSolutionLF = coarse_msfem_solution.local_discrete_function(coarse_entity);

    if (is_simplex_grid) {
      BOOST_ASSERT_MSG(localSolutions.size() == Dune::GridSelector::dimgrid,
                       "We should have dim local solutions per coarse element on triangular meshes!");

      MsFEMTraits::LocalGridDiscreteFunctionType::JacobianRangeType grad_coarse_msfem_on_entity;
      // We only need the gradient of the coarse scale part on the element, which is a constant.
      coarseSolutionLF.jacobian(coarse_entity.geometry().center(), grad_coarse_msfem_on_entity);

      // get the coarse gradient on T, multiply it with the local correctors and sum it up.
      for (int spaceDimension = 0; spaceDimension < Dune::GridSelector::dimgrid; ++spaceDimension) {
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
      for (std::size_t dof = 0; dof < coarseSolutionLF.vector().size(); ++dof) {
        localSolutions[dof]->vector() *= coarseSolutionLF.vector().get(dof);
        local_correction.vector() += localSolutions[dof]->vector();
      }

      // oversampling : restrict the local correctors to the element T
      // ie set all dofs not "covered" by the coarse cell to 0
      // also adds lg-prolongation of coarse_solution to local_correction
      const auto cut_overlay = DSC_CONFIG_GET("msfem.oversampling_layers", 0);
      for (auto& local_entity : DSC::viewRange(*localSolManager.space().grid_view())) {
        const auto& lg_points = localSolManager.space().lagrange_points(local_entity);
        const auto& reference_element = DSG::reference_element(coarse_entity);
        const auto& coarse_geometry = coarse_entity.geometry();
        auto entity_local_correction = local_correction.local_discrete_function(local_entity);

        for (const auto lg_i : DSC::valueRange(int(lg_points.size()))) {
          const auto local_point = lg_points[lg_i];
          const auto global_lg_point = local_entity.geometry().global(local_point);
          const auto local_coarse_point = coarse_geometry.local(global_lg_point);
          const auto coarse_value = coarseSolutionLF.evaluate(local_coarse_point) / double(lg_points.size());
          auto& vec = entity_local_correction.vector();
          if (cut_overlay) {
            const bool covered = reference_element.checkInside(local_coarse_point);
            vec.set(lg_i, covered ? vec.get(lg_i) : 0);
          }
          vec.add(lg_i, coarse_value);
        }
      }

      // add dirichlet corrector
      local_correction.vector() += localSolutions[coarseSolutionLF.vector().size() + 1]->vector();
      // substract neumann corrector
      local_correction.vector() -= localSolutions[coarseSolutionLF.vector().size()]->vector();

      if (DSC_CONFIG_GET("msfem.local_corrections_vtk_output", false)) {
        const std::string name = (boost::format("local_correction_%d_") % coarse_index).str();
        Dune::Multiscale::OutputParameters outputparam;
        outputparam.set_prefix(name);
        local_correction.visualize(outputparam.fullpath(local_correction.name()));
      }
    }
    localSolutions.clear();
  }

  msfem_solution = DSC::make_unique<LocalsolutionProxy>(std::move(local_corrections), coarse_space, subgrid_list);
}

void Elliptic_MsFEM_Solver::apply(const CommonTraits::SpaceType& coarse_space,
                                  std::unique_ptr<LocalsolutionProxy>& solution, LocalGridList& subgrid_list) const {
  DSC::Profiler::ScopedTiming st("msfem.Elliptic_MsFEM_Solver.apply");

  DiscreteFunctionType coarse_msfem_solution(coarse_space, "Coarse Part MsFEM Solution");
  coarse_msfem_solution.vector() *= 0;

  //! Solutions are kept in-memory via DiscreteFunctionIO::MemoryBackend by LocalsolutionManagers
  LocalProblemSolver(coarse_space, subgrid_list).solve_for_all_cells();

  CoarseScaleOperator elliptic_msfem_op(coarse_space, subgrid_list);
  elliptic_msfem_op.apply_inverse(coarse_msfem_solution);

  //! identify fine scale part of MsFEM solution (including the projection!)
  identify_fine_scale_part(subgrid_list, coarse_msfem_solution, solution);
}

} // namespace Multiscale {
} // namespace Dune {
