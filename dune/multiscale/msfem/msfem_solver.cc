#include <config.h>
#include <assert.h>
#include <boost/assert.hpp>
#include <dune/common/exceptions.hh>
#include <dune/common/timer.hh>
#include <dune/multiscale/common/righthandside_assembler.hh>
#include <dune/multiscale/msfem/coarse_scale_operator.hh>

#include <dune/stuff/fem/functions/checks.hh>
#include <dune/stuff/common/logging.hh>
#include <dune/stuff/common/parameter/configcontainer.hh>
#include <dune/stuff/common/profiler.hh>
#include <dune/stuff/discretefunction/projection/heterogenous.hh>
#include <dune/multiscale/msfem/localproblems/localgridsearch.hh>
#include <sstream>

#include <dune/multiscale/common/traits.hh>
#include <dune/multiscale/msfem/msfem_traits.hh>
#include <dune/gdt/operators/prolongations.hh>

#include "msfem_solver.hh"
#include "localsolution_proxy.hh"

namespace Dune {
namespace Multiscale {
namespace MsFEM {

void Elliptic_MsFEM_Solver::identify_fine_scale_part(LocalGridList& subgrid_list,
                                                     const DiscreteFunctionType& coarse_msfem_solution,
                                                     DiscreteFunctionType& fine_scale_part) const {
  DSC::Profiler::ScopedTiming st("msfem.idFine");
  fine_scale_part.vector() *= 0;
  const auto& coarse_space = coarse_msfem_solution.space();
  auto& coarse_indexset = coarse_space.grid_view()->grid().leafIndexSet();

  typedef DS::MsFEMProjection<LocalGridSearch> ProjectionType;
  typedef LocalGridSearch<typename MsFEMTraits::LocalGridViewType> SearchType;
  typedef LocalsolutionProxy<SearchType> ProxyType;
  typename ProxyType::CorrectionsMapType local_corrections;

  const bool is_simplex_grid = DSG::is_simplex_grid(coarse_space);

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
      BOOST_ASSERT_MSG(
          static_cast<long long>(localSolutions.size() - localSolManager.numBoundaryCorrectors()) ==
              static_cast<long long>(coarseSolutionLF.size()),
          "The current implementation relies on having thesame types of elements on coarse and fine level!");
      for (int dof = 0; dof < coarseSolutionLF.size(); ++dof) {
        localSolutions[dof]->vector() *= coarseSolutionLF.vector().get(dof);
        local_correction.vector() += localSolutions[dof]->vector();
      }

      // oversampling : restrict the local correctors to the element T
      // ie set all dofs not "covered" by the coarse cell to 0
      if (DSC_CONFIG_GET("msfem.oversampling_layers", 0)) {
        for (auto& local_entity : DSC::viewRange(*localSolManager.space().grid_view())) {
          DUNE_THROW(NotImplemented, "lg points in gdt?");
          const auto& lg_points = localSolManager.space().lagrange_points(local_entity);
          const auto& reference_element = DSG::reference_element(coarse_entity);
          const auto& coarse_geometry = coarse_entity.geometry();
          auto entity_local_correction = local_correction.local_discrete_function(local_entity);

          for (const auto lg_i : DSC::valueRange(int(lg_points.size()))) {
            const auto global_lg_point = local_entity.geometry().global(lg_points[lg_i]);
            const bool covered = reference_element.checkInside(coarse_geometry.local(global_lg_point));
            auto& vec = entity_local_correction.vector();
            vec.set(lg_i, covered ? vec.get(lg_i) : 0);
          }
        }
      }

      // add dirichlet corrector
      local_correction.vector() += localSolutions[coarseSolutionLF.vector().size() + 1]->vector();
      // substract neumann corrector
      local_correction.vector() -= localSolutions[coarseSolutionLF.vector().size()]->vector();

      if (DSC_CONFIG_GET("msfem.local_corrections_vtk_output", false)) {
        const auto elId = coarse_space.grid_view()->grid().globalIdSet().id(coarse_entity);
        const std::string name = (boost::format("local_correction_%d_") % elId).str();
        Dune::Multiscale::OutputParameters outputparam;
        outputparam.set_prefix(name);
        local_correction.visualize(outputparam.fullpath(local_correction.name()));
      }
    }
  }

  SearchType search(coarse_space, subgrid_list);
  ProxyType proxy(local_corrections, coarse_indexset, search);
  ProjectionType::project(proxy, fine_scale_part, search);
  BOOST_ASSERT_MSG(fine_scale_part.dofs_valid(), "Fine scale part DOFs need to be valid!");
  //backend storage no longer needed from here on
  DiscreteFunctionIO<MsFEMTraits::LocalGridDiscreteFunctionType>::clear();
}

void Elliptic_MsFEM_Solver::apply(const CommonTraits::DiscreteFunctionSpaceType& coarse_space,
                                  const CommonTraits::DiffusionType& diffusion_op,
                                  const CommonTraits::SourceType& f, DiscreteFunctionType& coarse_scale_part,
                                  DiscreteFunctionType& fine_scale_part, DiscreteFunctionType& solution) const {
  DSC::Profiler::ScopedTiming st("msfem.Elliptic_MsFEM_Solver.apply");
  BOOST_ASSERT_MSG(coarse_scale_part.dofs_valid(), "Coarse scale part DOFs need to be valid!");

  DiscreteFunctionType coarse_msfem_solution(coarse_space, "Coarse Part MsFEM Solution");
  coarse_msfem_solution.vector() *= 0;

  LocalGridList subgrid_list(coarse_space);
  //! Solutions are kept in-memory via DiscreteFunctionIO::MemoryBackend by LocalsolutionManagers
  LocalProblemSolver(coarse_space, subgrid_list, diffusion_op).solve_for_all_cells();

  DiscreteFunctionType msfem_rhs(coarse_space, "MsFEM right hand side");
  msfem_rhs.vector() *= 0;
  RightHandSideAssembler rhsAss;

//  rhsAss.assemble(coarse_space, f, subgrid_list, msfem_rhs);
  //! define the discrete (elliptic) operator that describes our problem
  CoarseScaleOperator elliptic_msfem_op(coarse_space, subgrid_list, diffusion_op);
  elliptic_msfem_op.apply_inverse(msfem_rhs, coarse_msfem_solution);

  solution.vector() *= 0;
//  Dune::Multiscale::copyDirichletValues(coarse_space, solution);
  DUNE_THROW(NotImplemented, "dirichlet ");

  //! identify fine scale part of MsFEM solution (including the projection!)
  identify_fine_scale_part(subgrid_list, coarse_msfem_solution, fine_scale_part);

  GDT::Operators::LagrangeProlongation< CommonTraits::GridViewType > projection(*coarse_msfem_solution.space().grid_view());
  projection.apply(coarse_msfem_solution, coarse_scale_part);
  // add coarse and fine scale part to solution
  solution.vector() += coarse_scale_part.vector();
  solution.vector() += fine_scale_part.vector();
} // solve_dirichlet_zero

} // namespace MsFEM {
} // namespace Multiscale {
} // namespace Dune {
