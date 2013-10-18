#include <config.h>
#include "msfem_solver.hh"

#include <unordered_set>

#include <dune/multiscale/common/righthandside_assembler.hh>
#include <dune/multiscale/msfem/localproblems/subgrid-list.hh>
#include <dune/multiscale/msfem/elliptic_msfem_matrix_assembler.hh>
#include <dune/multiscale/tools/misc/linear-lagrange-interpolation.hh>
#include <dune/multiscale/msfem/msfem_grid_specifier.hh>
#include <dune/multiscale/msfem/elliptic_msfem_matrix_assembler.hh>
#include <dune/multiscale/msfem/localproblems/localsolutionmanager.hh>
#include <dune/multiscale/fem/fem_traits.hh>
#include <dune/stuff/discretefunction/projection/heterogenous.hh>
#include <dune/stuff/common/profiler.hh>

#include <boost/assert.hpp>

namespace Dune {
namespace Multiscale {
namespace MsFEM {

Elliptic_MsFEM_Solver::Elliptic_MsFEM_Solver(const DiscreteFunctionSpace& discreteFunctionSpace)
  : discreteFunctionSpace_(discreteFunctionSpace) {}

void Elliptic_MsFEM_Solver::subgrid_to_hostrid_projection(const SubgridDiscreteFunctionType& sub_func,
                                                          DiscreteFunctionType& host_func) const {
  host_func.clear();

  const SubgridDiscreteFunctionSpaceType& subDiscreteFunctionSpace = sub_func.space();
  const SubGridType& subGrid = subDiscreteFunctionSpace.grid();

  typedef typename SubgridDiscreteFunctionSpaceType::IteratorType SubgridIterator;
  typedef typename SubgridIterator::Entity SubgridEntity;
  typedef typename SubgridDiscreteFunctionType::LocalFunctionType SubgridLocalFunction;

  const SubgridIterator sub_endit = subDiscreteFunctionSpace.end();
  for (SubgridIterator sub_it = subDiscreteFunctionSpace.begin(); sub_it != sub_endit; ++sub_it) {
    const SubgridEntity& sub_entity = *sub_it;

    const HostEntityPointer host_entity_pointer = subGrid.getHostEntity<0>(*sub_it);
    const HostEntity& host_entity = *host_entity_pointer;

    const SubgridLocalFunction sub_loc_value = sub_func.localFunction(sub_entity);
    auto host_loc_value = host_func.localFunction(host_entity);

    const auto numBaseFunctions = sub_loc_value.basisFunctionSet().size();
    for (unsigned int i = 0; i < numBaseFunctions; ++i) {
      host_loc_value[i] = sub_loc_value[i];
    }
  }
} // subgrid_to_hostrid_projection

void Elliptic_MsFEM_Solver::projectCoarseToFineScale(MacroMicroGridSpecifier& /*specifier*/,
                                                     const DiscreteFunctionType& coarse_msfem_solution,
                                                     DiscreteFunctionType& coarse_scale_part) const {

  DSC_LOG_INFO << "Indentifying coarse scale part of the MsFEM solution... ";

  coarse_scale_part.clear();
  Dune::Stuff::HeterogenousProjection<> projection;
  projection.project(coarse_msfem_solution, coarse_scale_part);
  DSC_LOG_INFO << " done." << std::endl;
}

void Elliptic_MsFEM_Solver::identify_fine_scale_part(MacroMicroGridSpecifier& specifier,
                                                     MsFEMTraits::SubGridListType& subgrid_list,
                                                     const DiscreteFunctionType& coarse_msfem_solution,
                                                     DiscreteFunctionType& fine_scale_part) const {
  fine_scale_part.clear();

  const GridPart& gridPart = discreteFunctionSpace_.gridPart();
  const HostGrid& grid = gridPart.grid();

  DiscreteFunctionSpace& coarse_space = specifier.coarseSpace();
  const HostGridLeafIndexSet& coarseGridLeafIndexSet = coarse_space.gridPart().grid().leafIndexSet();

  const int number_of_nodes = grid.size(2 /*codim*/);
  std::vector<std::vector<HostEntityPointer>> entities_sharing_same_node(number_of_nodes);

  DSC_LOG_INFO << "Indentifying fine scale part of the MsFEM solution... ";
  // traverse coarse space
  for (auto& coarseCell : coarse_space) {
    const auto coarseCellIndex = coarseGridLeafIndexSet.index(coarseCell);

    LocalSolutionManager localSolManager(coarseCell, subgrid_list, specifier);
    localSolManager.loadLocalSolutions();
    auto& localSolutions = localSolManager.getLocalSolutions();

    auto coarseSolutionLF = coarse_msfem_solution.localFunction(coarseCell);

    if ((specifier.getOversamplingStrategy() == 3) || specifier.simplexCoarseGrid()) {
      BOOST_ASSERT_MSG(localSolutions.size() == Dune::GridSelector::dimgrid,
                       "We should have dim local solutions per coarse element on triangular meshes!");

      JacobianRangeType grad_coarse_msfem_on_entity;
      // We only need the gradient of the coarse scale part on the element, which is a constant.
      coarseSolutionLF.jacobian(coarseCell.geometry().center(), grad_coarse_msfem_on_entity);

      // get the coarse gradient on T, multiply it with the local correctors and sum it up.
      for (int spaceDimension = 0; spaceDimension < Dune::GridSelector::dimgrid; ++spaceDimension) {
        *localSolutions[spaceDimension] *= grad_coarse_msfem_on_entity[0][spaceDimension];
        if (spaceDimension > 0)
          *localSolutions[0] += *localSolutions[spaceDimension];
      }
    } else {
      //! @warning At this point, we assume to have the same types of elements in the coarse and fine grid!
      BOOST_ASSERT_MSG(
          static_cast<long long>(localSolutions.size() - localSolManager.numBoundaryCorrectors()) ==
              static_cast<long long>(coarseSolutionLF.numDofs()),
          "The current implementation relies on having thesame types of elements on coarse and fine level!");
      for (int dof = 0; dof < coarseSolutionLF.numDofs(); ++dof) {
        *localSolutions[dof] *= coarseSolutionLF[dof];
        if (dof > 0)
          *localSolutions[0] += *localSolutions[dof];
      }

      // add dirichlet corrector
      DiscreteFunctionType boundaryCorrector("Boundary Corrector", discreteFunctionSpace_);
      boundaryCorrector.clear();
      subgrid_to_hostrid_projection(*localSolutions[coarseSolutionLF.numDofs() + 1], boundaryCorrector);
      fine_scale_part += boundaryCorrector;

      // substract neumann corrector
      // boundaryCorrector.clear();
      // subgrid_to_hostrid_projection(*localSolutions[coarseSolutionLF.numDofs()], boundaryCorrector);
      // fine_scale_part -= boundaryCorrector;
    }

    // oversampling strategy 3: just sum up the local correctors:
    if ((specifier.getOversamplingStrategy() == 3)) {
      DiscreteFunctionType correction_on_U_T("correction_on_U_T", discreteFunctionSpace_);
      subgrid_to_hostrid_projection(*localSolutions[0], correction_on_U_T);
      fine_scale_part += correction_on_U_T;
    }

    // oversampling strategy 1 or 2: restrict the local correctors to the element T, sum them up and apply a conforming
    // projection:
    if ((specifier.getOversamplingStrategy() == 1) || (specifier.getOversamplingStrategy() == 2)) {
      BOOST_ASSERT_MSG(localSolManager.getLocalDiscreteFunctionSpace().gridPart().grid().maxLevel() ==
                           discreteFunctionSpace_.gridPart().grid().maxLevel(),
                       "Error: MaxLevel of SubGrid not identical to MaxLevel of FineGrid.");

      const auto& nodeToEntityMap = subgrid_list.getNodeEntityMap();

      for (auto& subgridEntity : localSolManager.getLocalDiscreteFunctionSpace()) {
        //! MARK actual subgrid usage
        const auto fine_host_entity_pointer = localSolManager.getSubGridPart().grid().getHostEntity<0>(subgridEntity);
        const auto& fine_host_entity = *fine_host_entity_pointer;

        const auto hostFatherIndex = subgrid_list.getEnclosingMacroCellIndex(fine_host_entity_pointer);
        if (hostFatherIndex == coarseCellIndex) {
          const auto sub_loc_value = localSolutions[0]->localFunction(subgridEntity);

          assert(localSolutions.size() == coarseSolutionLF.numDofs() + localSolManager.numBoundaryCorrectors());
          auto host_loc_value = fine_scale_part.localFunction(fine_host_entity);

          const auto number_of_nodes_entity = subgridEntity.count<HostGrid::dimension>();

          for (auto i : DSC::valueRange(number_of_nodes_entity)) {
            const auto node = fine_host_entity.subEntity<HostGrid::dimension>(i);
            const auto global_index_node = gridPart.grid().leafIndexSet().index(*node);

            // count the number of different coarse-grid-entities that share the above node
            std::unordered_set<SubGridListType::IdType> coarse_entities;
            const auto numEntitiesSharingNode = nodeToEntityMap[global_index_node].size();
            for (size_t j = 0; j < numEntitiesSharingNode; ++j) {
              // get the id of the macro element enclosing the current element
              const auto innerId = subgrid_list.getEnclosingMacroCellId(nodeToEntityMap[global_index_node][j]);
              // the following will only add the entity index if it is not yet present
              coarse_entities.insert(innerId);
            }
            host_loc_value[i] += (sub_loc_value[i] / coarse_entities.size());
          }
        }
      }
    }
  }
  DSC_LOG_INFO << " done." << std::endl;
}

void Elliptic_MsFEM_Solver::solve_dirichlet_zero(
    const CommonTraits::DiffusionType& diffusion_op, const CommonTraits::FirstSourceType& f,
    // number of layers per coarse grid entity T:  U(T) is created by enrichting T with
    // n(T)-layers.
    MacroMicroGridSpecifier& specifier, MsFEMTraits::SubGridListType& subgrid_list,
    DiscreteFunctionType& coarse_scale_part, DiscreteFunctionType& fine_scale_part,
    DiscreteFunctionType& solution) const {
  DSC::Profiler::ScopedTiming st("msfem.Elliptic_MsFEM_Solver.solve_dirichlet_zero");

  DiscreteFunctionSpace& coarse_space = specifier.coarseSpace();

  DiscreteFunctionType coarse_msfem_solution("Coarse Part MsFEM Solution", coarse_space);
  coarse_msfem_solution.clear();

  //! define the right hand side assembler tool
  // (for linear and non-linear elliptic and parabolic problems, for sources f and/or G )
  typedef RightHandSideAssembler RhsAssembler;

  //! define the discrete (elliptic) operator that describes our problem
  // discrete elliptic MsFEM operator (corresponds with MsFEM Matrix)
  // ( effect of the discretized differential operator on a certain discrete function )
  // This will assemble and solve the local problems
  const DiscreteEllipticMsFEMOperator elliptic_msfem_op(specifier, coarse_space, subgrid_list, diffusion_op);
  // discrete elliptic operator (corresponds with FEM Matrix)

  //! (stiffness) matrix
  MsLinearOperatorTypeType msfem_matrix("MsFEM stiffness matrix", coarse_space, coarse_space);

  //! right hand side vector
  // right hand side for the finite element method:
  DiscreteFunctionType msfem_rhs("MsFEM right hand side", coarse_space);
  msfem_rhs.clear();

  DSC_LOG_INFO << std::endl << "Solving MsFEM problem." << std::endl
               << "Solving linear problem with MsFEM and maximum coarse grid level "
               << coarse_space.gridPart().grid().maxLevel() << "." << std::endl
               << "------------------------------------------------------------------------------" << std::endl;

  // to assemble the computational time
  Dune::Timer assembleTimer;

  // assemble the MsFEM stiffness matrix
  elliptic_msfem_op.assemble_matrix(msfem_matrix);

  DSC_LOG_INFO << "Time to assemble MsFEM stiffness matrix: " << assembleTimer.elapsed() << "s" << std::endl;

  // assemble right hand side
  if (DSC_CONFIG_GET("msfem.petrov_galerkin", 1)) {
    RhsAssembler::assemble(f, msfem_rhs);
  } else {
    RhsAssembler::assemble_for_MsFEM_symmetric(f, specifier, subgrid_list, msfem_rhs);
  }
  msfem_rhs.communicate();
  BOOST_ASSERT_MSG(msfem_rhs.dofsValid(), "Coarse scale RHS DOFs need to be valid!");

  const InverseOperatorType msfem_biCGStab(msfem_matrix, 1e-8, 1e-8, 2000, true, "bcgs",
                                           DSC_CONFIG_GET("preconditioner_type", std::string("sor")));
  msfem_biCGStab(msfem_rhs, coarse_msfem_solution);
  DSC_LOG_INFO << "---------------------------------------------------------------------------------" << std::endl;
  DSC_LOG_INFO << "MsFEM problem solved in " << assembleTimer.elapsed() << "s." << std::endl << std::endl << std::endl;

  if (!coarse_msfem_solution.dofsValid())
    DUNE_THROW(InvalidStateException, "Degrees of freedom of coarse solution are not valid!");

  // get the dirichlet values
  solution.clear();
  Dune::Multiscale::projectDirichletValues(coarse_space, solution);

  //! copy coarse scale part of MsFEM solution into a function defined on the fine grid
  projectCoarseToFineScale(specifier, coarse_msfem_solution, coarse_scale_part);

  //! identify fine scale part of MsFEM solution (including the projection!)
  identify_fine_scale_part(specifier, subgrid_list, coarse_msfem_solution, fine_scale_part);
  {
    DSC::Profiler::ScopedTiming commFSTimer("msfem.Elliptic_MsFEM_Solver.solve_dirichlet_zero.comm_fine_scale_part");
    fine_scale_part.communicate();
  }

  BOOST_ASSERT_MSG(coarse_scale_part.dofsValid(), "Coarse scale part DOFs need to be valid!");
  BOOST_ASSERT_MSG(fine_scale_part.dofsValid(), "Fine scale part DOFs need to be valid!");

  // add coarse and fine scale part to solution
  solution += coarse_scale_part;
  solution += fine_scale_part;

} // solve_dirichlet_zero

} // namespace MsFEM {
} // namespace Multiscale {
} // namespace Dune {
