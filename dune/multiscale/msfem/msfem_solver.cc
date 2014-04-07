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

#include "dune/multiscale/common/dirichletconstraints.hh"
#include "dune/multiscale/common/traits.hh"
#include "dune/multiscale/msfem/msfem_traits.hh"
#include <dune/multiscale/common/output_traits.hh>
#include "msfem_solver.hh"
#include "localsolution_proxy.hh"

namespace Dune {
namespace Multiscale {
namespace MsFEM {

struct LocalOutputTraits {
  //! --------- typedefs and classes for data output -----------------------------------------
  typedef std::tuple<const MsFEMTraits::LocalGridDiscreteFunctionType*> IOTupleType;
  typedef Dune::Fem::DataWriter<MsFEMTraits::LocalGridType, IOTupleType> DataOutputType;
};

void Elliptic_MsFEM_Solver::identify_fine_scale_part(LocalGridList& subgrid_list,
                                                     const DiscreteFunctionType& coarse_msfem_solution,
                                                     DiscreteFunctionType& fine_scale_part) const {
  fine_scale_part.clear();
  const DiscreteFunctionSpace& coarse_space = coarse_msfem_solution.space();
  auto& coarse_indexset = coarse_space.gridPart().grid().leafIndexSet();

  typedef DS::HeterogenousProjection<LocalGridSearch> ProjectionType;
  typedef LocalGridSearch<typename MsFEMTraits::LocalGridType::LeafGridView> SearchType;
  typedef LocalsolutionProxy<SearchType> ProxyType;
  typename ProxyType::CorrectionsMapType local_corrections;

  for (auto& coarse_entity : coarse_space) {
    LocalSolutionManager localSolManager(coarse_space, coarse_entity, subgrid_list);
    localSolManager.load();
    auto& localSolutions = localSolManager.getLocalSolutions();
    auto coarse_index = coarse_indexset.index(coarse_entity);
    local_corrections[coarse_index] =
        DSC::make_unique<MsFEMTraits::LocalGridDiscreteFunctionType>("", localSolManager.space());

    auto& local_correction = *local_corrections[coarse_index];
    local_correction.clear();
    auto coarseSolutionLF = coarse_msfem_solution.localFunction(coarse_entity);

    if (DSG::is_simplex_grid(coarse_space)) {
      BOOST_ASSERT_MSG(localSolutions.size() == Dune::GridSelector::dimgrid,
                       "We should have dim local solutions per coarse element on triangular meshes!");

      JacobianRangeType grad_coarse_msfem_on_entity;
      // We only need the gradient of the coarse scale part on the element, which is a constant.
      coarseSolutionLF.jacobian(coarse_entity.geometry().center(), grad_coarse_msfem_on_entity);

      // get the coarse gradient on T, multiply it with the local correctors and sum it up.
      for (int spaceDimension = 0; spaceDimension < Dune::GridSelector::dimgrid; ++spaceDimension) {
        *localSolutions[spaceDimension] *= grad_coarse_msfem_on_entity[0][spaceDimension];
        local_correction += *localSolutions[spaceDimension];
      }
      BOOST_ASSERT_MSG(false, "no adding of the boundary correctors??");
    } else {
      //! @warning At this point, we assume to have the same types of elements in the coarse and fine grid!
      BOOST_ASSERT_MSG(
          static_cast<long long>(localSolutions.size() - localSolManager.numBoundaryCorrectors()) ==
              static_cast<long long>(coarseSolutionLF.numDofs()),
          "The current implementation relies on having thesame types of elements on coarse and fine level!");
      for (int dof = 0; dof < coarseSolutionLF.numDofs(); ++dof) {
        *localSolutions[dof] *= coarseSolutionLF[dof];
        local_correction += *localSolutions[dof];
      }

      // oversampling : restrict the local correctors to the element T
      // ie set all dofs not "covered" by the coarse cell to 0
      if (DSC_CONFIG_GET("msfem.oversampling_layers", 0)) {
        for (auto& local_entity : localSolManager.space()) {
          const auto& lg_points = localSolManager.space().lagrangePointSet(local_entity);
          const auto& reference_element = DSG::reference_element(coarse_entity);
          const auto& coarse_geometry = coarse_entity.geometry();
          auto entity_local_correction = local_correction.localFunction(local_entity);

          for (const auto lg_i : DSC::valueRange(int(lg_points.size()))) {
            const auto global_lg_point = local_entity.geometry().global(lg_points.point(lg_i));
            const bool covered = reference_element.checkInside(coarse_geometry.local(global_lg_point));
            entity_local_correction[lg_i] = covered ? entity_local_correction[lg_i] : RangeType::value_type(0);
          }
        }
      }

      // add dirichlet corrector
      local_correction += *localSolutions[coarseSolutionLF.numDofs() + 1];
      // substract neumann corrector
      local_correction -= *localSolutions[coarseSolutionLF.numDofs()];

      if (DSC_CONFIG_GET("msfem.local_corrections_vtk_output", false)) {
        LocalOutputTraits::IOTupleType coarse_grid_series(&local_correction);
        const auto elId = coarse_space.gridPart().grid().globalIdSet().id(coarse_entity);
        const std::string name = (boost::format("local_correction_%d_") % elId).str();
        Dune::Multiscale::OutputParameters outputparam;
        outputparam.set_prefix(name);
        auto& grid = localSolManager.space().gridPart().grid();
        LocalOutputTraits::DataOutputType(grid, coarse_grid_series, outputparam).writeData(1.0 /*dummy*/, name);
      }
    }
  }

  SearchType search(coarse_space, subgrid_list);
  ProxyType proxy(local_corrections, coarse_indexset, search);
  ProjectionType::project(proxy, fine_scale_part, search);
  BOOST_ASSERT_MSG(fine_scale_part.dofsValid(), "Fine scale part DOFs need to be valid!");
  //backend storage no longer needed from here on
  DiscreteFunctionIO<MsFEMTraits::LocalGridDiscreteFunctionType>::clear();
}

void Elliptic_MsFEM_Solver::apply(const CommonTraits::DiscreteFunctionSpaceType& coarse_space,
                                  const CommonTraits::DiffusionType& diffusion_op,
                                  const CommonTraits::FirstSourceType& f, DiscreteFunctionType& coarse_scale_part,
                                  DiscreteFunctionType& fine_scale_part, DiscreteFunctionType& solution) const {
  if (DSC_CONFIG_GET("msfem.petrov_galerkin", 1))
    DSC_LOG_ERROR << "MsFEM does not work with Petrov-Galerkin at the moment!\n";

  DSC::Profiler::ScopedTiming st("msfem.Elliptic_MsFEM_Solver.apply");
  BOOST_ASSERT_MSG(coarse_scale_part.dofsValid(), "Coarse scale part DOFs need to be valid!");

  DiscreteFunctionType coarse_msfem_solution("Coarse Part MsFEM Solution", coarse_space);
  coarse_msfem_solution.clear();

  LocalGridList subgrid_list(coarse_space);
  //! Solutions are kept in-memory via DiscreteFunctionIO::MemoryBackend by LocalsolutionManagers
  LocalProblemSolver(coarse_space, subgrid_list, diffusion_op).solve_for_all_cells();

  DiscreteFunctionType msfem_rhs("MsFEM right hand side", coarse_space);
  msfem_rhs.clear();
  RightHandSideAssembler::assemble_msfem(coarse_space, f, subgrid_list, msfem_rhs);
  //! define the discrete (elliptic) operator that describes our problem
  CoarseScaleOperator elliptic_msfem_op(coarse_space, subgrid_list, diffusion_op);
  elliptic_msfem_op.apply_inverse(msfem_rhs, coarse_msfem_solution);

  solution.clear();
  Dune::Multiscale::copyDirichletValues(coarse_space, solution);

  //! identify fine scale part of MsFEM solution (including the projection!)
  identify_fine_scale_part(subgrid_list, coarse_msfem_solution, fine_scale_part);

  fine_scale_part.communicate();

  DS::HeterogenousProjection<>::project(coarse_msfem_solution, coarse_scale_part);
  // add coarse and fine scale part to solution
  solution += coarse_scale_part;
  solution += fine_scale_part;

  // seperate the msfem output from other output
  std::cout << std::endl << std::endl;
} // solve_dirichlet_zero

} // namespace MsFEM {
} // namespace Multiscale {
} // namespace Dune {
