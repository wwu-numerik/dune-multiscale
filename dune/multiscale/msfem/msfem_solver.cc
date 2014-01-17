#include <config.h>
#include <assert.h>
#include <boost/assert.hpp>
#include <dune/common/exceptions.hh>
#include <dune/common/timer.hh>
#include <dune/multiscale/common/righthandside_assembler.hh>
#include <dune/multiscale/msfem/elliptic_msfem_matrix_assembler.hh>
#include <dune/multiscale/msfem/msfem_grid_specifier.hh>
#include <dune/stuff/common/logging.hh>
#include <dune/stuff/common/parameter/configcontainer.hh>
#include <dune/stuff/common/profiler.hh>
#include <dune/stuff/discretefunction/projection/heterogenous.hh>
#include <sstream>

#include "dune/multiscale/common/dirichletconstraints.hh"
#include "dune/multiscale/common/traits.hh"
#include "dune/multiscale/msfem/msfem_traits.hh"
#include "msfem_solver.hh"

namespace Dune {
namespace Multiscale {
namespace MsFEM {

Elliptic_MsFEM_Solver::Elliptic_MsFEM_Solver(const DiscreteFunctionSpace& discreteFunctionSpace)
  : discreteFunctionSpace_(discreteFunctionSpace) {}

void Elliptic_MsFEM_Solver::subgrid_to_hostrid_projection(const SubgridDiscreteFunctionType& sub_func,
                                                          DiscreteFunctionType& host_func) const {
  host_func.clear();
  Stuff::HeterogenousProjection<>::project(sub_func, host_func);
} // subgrid_to_hostrid_projection

void Elliptic_MsFEM_Solver::identify_fine_scale_part(MacroMicroGridSpecifier& specifier,
                                                     MsFEMTraits::LocalGridListType& subgrid_list,
                                                     const DiscreteFunctionType& coarse_msfem_solution,
                                                     DiscreteFunctionType& fine_scale_part) const {
  fine_scale_part.clear();
  DiscreteFunctionSpace& coarse_space = specifier.coarseSpace();

  DSC_LOG_INFO << "Indentifying fine scale part of the MsFEM solution... ";
  DiscreteFunctionType fine_correction("Boundary Corrector", discreteFunctionSpace_);

  // traverse coarse space
  for (auto& coarseCell : coarse_space) {
    LocalSolutionManager localSolManager(coarseCell, subgrid_list, specifier);
    localSolManager.loadLocalSolutions();
    auto& localSolutions = localSolManager.getLocalSolutions();

    auto coarseSolutionLF = coarse_msfem_solution.localFunction(coarseCell);
    DMM::MsFEMTraits::LocalGridDiscreteFunctionType local_correction("", localSolManager.space());
    local_correction.clear();

    if ((DSC_CONFIG_GET("msfem.oversampling_strategy", 1) == 3) || specifier.simplexCoarseGrid()) {
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
      local_correction += *localSolutions[coarseSolutionLF.numDofs() + 1];
      // substract neumann corrector
      local_correction -= *localSolutions[coarseSolutionLF.numDofs() + 1];
    }

    // oversampling strategy 3: just sum up the local correctors:
    if ((DSC_CONFIG_GET("msfem.oversampling_strategy", 1) == 3)) {
      local_correction += *localSolutions[0];
    }

    // oversampling strategy 1 or 2: restrict the local correctors to the element T, sum them up and apply a conforming
    // projection:
    if ((DSC_CONFIG_GET("msfem.oversampling_strategy", 1) == 1) || (DSC_CONFIG_GET("msfem.oversampling_strategy", 1) == 2)) {

//      DUNE_THROW(NotImplemented, "pretty sure this is bs. there's no sum of local correctors. restriction also no longer works");
//      for (auto& local_entity : localSolManager.space()) {
//        if (subgrid_list.covers(coarseCell, local_entity)) {
//          const auto sub_loc_value = localSolutions[0]->localFunction(local_entity);

//          assert(localSolutions.size() == coarseSolutionLF.numDofs() + localSolManager.numBoundaryCorrectors());
//          auto host_loc_value = fine_scale_part.localFunction(local_entity);

//          const auto number_of_nodes_entity = local_entity.count<LocalGrid::dimension>();

//          for (auto i : DSC::valueRange(number_of_nodes_entity)) {
//            const auto node = local_entity.subEntity<LocalGrid::dimension>(i);
//            const auto global_index_node = gridPart.indexSet().index(*node);

//            // devide the value by the number of fine elements sharing the node (will be
//            // added numEntitiesSharingNode times)
//            const auto numEntitiesSharingNode = nodeToEntityMap[global_index_node].size();
//            host_loc_value[i] += (sub_loc_value[i] / numEntitiesSharingNode);
//          }
//        }
//      }
    }

    Stuff::HeterogenousProjection<>::project(local_correction, fine_correction);
    fine_scale_part += fine_correction;
  }
  DSC_LOG_INFO << " done." << std::endl;
}

void Elliptic_MsFEM_Solver::solve_dirichlet_zero(
    const CommonTraits::DiffusionType& diffusion_op, const CommonTraits::FirstSourceType& f,
    // number of layers per coarse grid entity T:  U(T) is created by enrichting T with
    // n(T)-layers.
    MacroMicroGridSpecifier& specifier, MsFEMTraits::LocalGridListType& subgrid_list,
    DiscreteFunctionType& coarse_scale_part, DiscreteFunctionType& fine_scale_part,
    DiscreteFunctionType& solution) const {
  DSC::Profiler::ScopedTiming st("msfem.Elliptic_MsFEM_Solver.solve_dirichlet_zero");

  DiscreteFunctionSpace& coarse_space = specifier.coarseSpace();

  DiscreteFunctionType coarse_msfem_solution("Coarse Part MsFEM Solution", coarse_space);
  coarse_msfem_solution.clear();

  //! define the right hand side assembler tool
  // (for linear and non-linear elliptic and parabolic problems, for sources f and/or G )
  typedef RightHandSideAssembler RhsAssembler;

  // Assemble and solve the local problems. Timing is done in assembleAndSolveAll-method
  MsFEMLocalProblemSolver localProblemSolver(specifier, subgrid_list, diffusion_op);
  localProblemSolver.assembleAndSolveAll();

  //! define the discrete (elliptic) operator that describes our problem
  // discrete elliptic MsFEM operator (corresponds with MsFEM Matrix)
  // ( effect of the discretized differential operator on a certain discrete function )
  const DiscreteEllipticMsFEMOperator elliptic_msfem_op(specifier, coarse_space, subgrid_list, diffusion_op);

  //! (stiffness) matrix
  MsLinearOperatorTypeType msfem_matrix("MsFEM stiffness matrix", coarse_space, coarse_space);

  //! right hand side vector
  // right hand side for the finite element method:
  DiscreteFunctionType msfem_rhs("MsFEM right hand side", coarse_space);
  msfem_rhs.clear();

  // to assemble the computational time
  DSC_PROFILER.startTiming("msfem.assembleMatrix");
  // assemble the MsFEM stiffness matrix
  elliptic_msfem_op.assemble_matrix(msfem_matrix);
  DSC_LOG_INFO << "Time to assemble MsFEM stiffness matrix: " << DSC_PROFILER.stopTiming("msfem.assembleMatrix")
               << "ms" << std::endl;

  // assemble right hand side
  DSC_PROFILER.startTiming("msfem.assembleRHS");
  if (DSC_CONFIG_GET("msfem.petrov_galerkin", 1))
    DSC_LOG_ERROR << "MsFEM does not work with Petrov-Galerkin at the moment!\n";

  RhsAssembler::assemble_for_MsFEM_symmetric(f, specifier, subgrid_list, msfem_rhs);

  msfem_rhs.communicate();
  DSC_LOG_INFO << "Time to assemble and communicate MsFEM rhs: " << DSC_PROFILER.stopTiming("msfem.assembleRHS")
               << "ms" << std::endl;

  BOOST_ASSERT_MSG(msfem_rhs.dofsValid(), "Coarse scale RHS DOFs need to be valid!");

  DSC_PROFILER.startTiming("msfem.solveCoarse");
  const InverseOperatorType msfem_biCGStab(msfem_matrix, 1e-8, 1e-8,
                                           DSC_CONFIG_GET("msfem.solver.iterations", msfem_rhs.size()),
                                           DSC_CONFIG_GET("msfem.solver.verbose", false), "bcgs",
                                           DSC_CONFIG_GET("msfem.solver.preconditioner_type", std::string("sor")));
  msfem_biCGStab(msfem_rhs, coarse_msfem_solution);
  DSC_LOG_INFO << "Time to solve coarse MsFEM problem: " << DSC_PROFILER.stopTiming("msfem.solveCoarse")
               << "ms." << std::endl;

  if (!coarse_msfem_solution.dofsValid())
    DUNE_THROW(InvalidStateException, "Degrees of freedom of coarse solution are not valid!");

  // get the dirichlet values
  solution.clear();
  Dune::Multiscale::copyDirichletValues(coarse_space, solution);

  //! identify fine scale part of MsFEM solution (including the projection!)
  identify_fine_scale_part(specifier, subgrid_list, coarse_msfem_solution, fine_scale_part);
  {
    DSC::Profiler::ScopedTiming commFSTimer("msfem.Elliptic_MsFEM_Solver.solve_dirichlet_zero.comm_fine_scale_part");
    fine_scale_part.communicate();
  }

  BOOST_ASSERT_MSG(coarse_scale_part.dofsValid(), "Coarse scale part DOFs need to be valid!");
  BOOST_ASSERT_MSG(fine_scale_part.dofsValid(), "Fine scale part DOFs need to be valid!");

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
