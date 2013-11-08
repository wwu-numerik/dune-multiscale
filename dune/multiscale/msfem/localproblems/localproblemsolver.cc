#include <config.h>
#include <assert.h>
#include <boost/concept/usage.hpp>
#include <boost/format.hpp>
#include <dune/common/exceptions.hh>
#include <dune/fem/io/file/dataoutput.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/istl/scalarproducts.hh>
#include <dune/istl/solvers.hh>
#include <dune/multiscale/common/dirichletconstraints.hh>
#include <dune/multiscale/msfem/localproblems/localoperator.hh>
#include <dune/multiscale/problems/selector.hh>
#include <dune/multiscale/tools/misc/uzawa.hh>
#include <dune/multiscale/tools/discretefunctionwriter.hh>
#include <dune/multiscale/msfem/localproblems/localsolutionmanager.hh>
#include <dune/stuff/common/logging.hh>
#include <dune/stuff/common/math.hh>
#include <dune/stuff/common/memory.hh>
#include <dune/stuff/common/parameter/configcontainer.hh>
#include <dune/stuff/common/profiler.hh>
#include <dune/stuff/common/ranges.hh>
#include <dune/stuff/fem/localmatrix_proxy.hh>
#include <iterator>
#include <memory>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include "dune/multiscale/msfem/localproblems/subgrid-list.hh"
#include "dune/multiscale/tools/misc/outputparameter.hh"
#include "localproblemsolver.hh"

namespace Dune {
namespace Multiscale {
namespace MsFEM {

LocalProblemDataOutputParameters::LocalProblemDataOutputParameters()
  : OutputParameters(DSC_CONFIG_GET("global.datadir", "data") + "/local_problems/") {}

std::unique_ptr<MsFEMLocalProblemSolver::InverseLocProbLinearOperatorTypeType>
MsFEMLocalProblemSolver::make_inverse_operator(MsFEMLocalProblemSolver::LocProbLinearOperatorTypeType& problem_matrix) {
  const auto solver =
      Dune::Multiscale::Problem::getModelData()->symmetricDiffusion() ? std::string("cg") : std::string("bcgs");
  return DSC::make_unique<InverseLocProbLinearOperatorTypeType>(
      problem_matrix, 1e-8, 1e-8, 20000, DSC_CONFIG_GET("localproblemsolver_verbose", false), solver,
      DSC_CONFIG_GET("preconditioner_type", std::string("sor")), 1);
}

MsFEMLocalProblemSolver::MsFEMLocalProblemSolver(const HostDiscreteFunctionSpaceType& hostDiscreteFunctionSpace,
                                                 const MacroMicroGridSpecifierType& specifier,
                                                 SubGridList& subgrid_list,
                                                 const DiffusionOperatorType& diffusion_operator)
  : hostDiscreteFunctionSpace_(hostDiscreteFunctionSpace)
  , diffusion_(diffusion_operator)
  , specifier_(specifier)
  , subgrid_list_(subgrid_list)
  , ids_relevant_basis_functions_for_subgrid_(nullptr)
  , inverse_of_L1_norm_coarse_basis_funcs_(nullptr)
  , coarse_basis_(nullptr)
  , global_id_to_internal_id_(nullptr)
  , neumann_bc_(nullptr)
  , dirichlet_extension_(nullptr)
  , threadIterators_(specifier_.coarseSpace().gridPart()) {}

MsFEMLocalProblemSolver::MsFEMLocalProblemSolver(
    const HostDiscreteFunctionSpaceType& hostDiscreteFunctionSpace, const MacroMicroGridSpecifierType& specifier,
    SubGridList& subgrid_list, std::vector<std::vector<std::size_t>>& ids_basis_functions_in_subgrid,
    std::vector<double>& inverse_of_L1_norm_coarse_basis_funcs, const DiffusionOperatorType& diffusion_operator,
    const CoarseBasisFunctionListType& coarse_basis, const std::map<std::size_t, std::size_t>& global_id_to_internal_id,
    const NeumannBoundaryType& neumann_bc, const HostDiscreteFunctionType& dirichlet_extension)
  : hostDiscreteFunctionSpace_(hostDiscreteFunctionSpace)
  , diffusion_(diffusion_operator)
  , specifier_(specifier)
  , subgrid_list_(subgrid_list)
  , ids_relevant_basis_functions_for_subgrid_(&ids_basis_functions_in_subgrid)
  ,
  // ids of the coarse grid basis functions in the interior of the subgrid
  inverse_of_L1_norm_coarse_basis_funcs_(&inverse_of_L1_norm_coarse_basis_funcs)
  , coarse_basis_(&coarse_basis)
  , global_id_to_internal_id_(&global_id_to_internal_id)
  , neumann_bc_(&neumann_bc)
  , dirichlet_extension_(&dirichlet_extension)
  , threadIterators_(specifier_.coarseSpace().gridPart()) {}

//! ----------- method: solve the local MsFEM problem ------------------------------------------
/** Solve all local MsFEM problems for one coarse entity at once.
*
*
*/
void MsFEMLocalProblemSolver::solveAllLocalProblems(const CoarseEntityType& coarseCell,
                                                    SubDiscreteFunctionVectorType& allLocalSolutions) const {
  assert(allLocalSolutions.size() > 0);

  const bool hasBoundary = coarseCell.hasBoundaryIntersections();
  const auto numBoundaryCorrectors = specifier_.simplexCoarseGrid() ? 1u : 2u;
  const auto numInnerCorrectors = allLocalSolutions.size() - numBoundaryCorrectors;

  // clear return argument
  for (auto& localSol : allLocalSolutions)
    localSol->clear();

  const auto& subDiscreteFunctionSpace = allLocalSolutions[0]->space();

  //! the matrix in our linear system of equations
  // in the non-linear case, it is the matrix for each iteration step
  LocProbLinearOperatorTypeType locProbSysMatrix("Local Problem System Matrix", subDiscreteFunctionSpace,
                                                 subDiscreteFunctionSpace);

  //! define the discrete (elliptic) local MsFEM problem operator
  // ( effect of the discretized differential operator on a certain discrete function )
  LocalProblemOperator localProblemOperator(subDiscreteFunctionSpace, diffusion_);

  // right hand side vector of the algebraic local MsFEM problem
  SubDiscreteFunctionVectorType allLocalRHS(allLocalSolutions.size());
  for (auto& it : allLocalRHS)
    it = DSC::make_unique<SubDiscreteFunctionType>("rhs of local MsFEM problem", subDiscreteFunctionSpace);

  switch (specifier_.getOversamplingStrategy()) {
    case 1:
      localProblemOperator.assemble_matrix(locProbSysMatrix);
      localProblemOperator.assembleAllLocalRHS(coarseCell, specifier_, allLocalRHS);
      break;
    default:
      DUNE_THROW(Fem::ParameterInvalid, "Oversampling Strategy must be 1 at the moment");
  }

  // set dirichlet dofs to zero
  Stuff::GridboundaryAllDirichlet<SubGridType::LeafGridView::Intersection> boundaryInfo;
  DirichletConstraints<SubDiscreteFunctionSpaceType> constraints(boundaryInfo, subDiscreteFunctionSpace);
  constraints.applyToOperator(locProbSysMatrix);

  for (auto& rhsIt : allLocalRHS) {
    constraints.setValue(0.0, *rhsIt);
  }
  for (auto i : DSC::valueRange(allLocalSolutions.size())) {
    if (!allLocalRHS[i]->dofsValid())
      DUNE_THROW(Dune::InvalidStateException, "Local MsFEM Problem RHS invalid.");

    // is the right hand side of the local MsFEM problem equal to zero or almost identical to zero?
    // if yes, the solution of the local MsFEM problem is also identical to zero. The solver is getting a problem with
    // this situation, which is why we do not solve local msfem problems for zero-right-hand-side, since we already know
    // the result.
    if (localProblemOperator.normRHS(*allLocalRHS[i]) < 1e-30) {
      allLocalRHS[i]->clear();
      DSC_LOG_ERROR << "Local MsFEM problem with solution zero." << std::endl;
      continue;
    }

    // don't solve local problems for boundary correctors if coarse cell has no boundary intersections
    if (i >= numInnerCorrectors && !hasBoundary) {
      allLocalRHS[i]->clear();
      DSC_LOG_INFO << "Zero-Boundary corrector." << std::endl;
      continue;
    }

    const auto localProblemSolver = make_inverse_operator(locProbSysMatrix);
    localProblemSolver->apply(*allLocalRHS[i], *allLocalSolutions[i]);

    if (!(allLocalSolutions[i]->dofsValid()))
      DUNE_THROW(Dune::InvalidStateException, "Current solution of the local msfem problem invalid!");
  }

  return;
}

//! ----------- method: solve the local MsFEM problem ------------------------------------------

void MsFEMLocalProblemSolver::solvelocalproblem(JacobianRangeType& e, SubDiscreteFunctionType& local_problem_solution,
                                                const int coarse_index /*= -1*/) const {
  // set solution equal to zero:
  local_problem_solution.clear();

  const auto& subDiscreteFunctionSpace = local_problem_solution.space();

  //! the matrix in our linear system of equations
  // in the non-linear case, it is the matrix for each iteration step
  LocProbLinearOperatorTypeType locprob_system_matrix("Local Problem System Matrix", subDiscreteFunctionSpace,
                                                      subDiscreteFunctionSpace);

  //! define the discrete (elliptic) local MsFEM problem operator
  // ( effect of the discretized differential operator on a certain discrete function )
  LocalProblemOperator local_problem_op(subDiscreteFunctionSpace, diffusion_);

  //! right hand side vector of the algebraic local MsFEM problem
  SubDiscreteFunctionType local_problem_rhs("rhs of local MsFEM problem", subDiscreteFunctionSpace);
  local_problem_rhs.clear();
  switch (specifier_.getOversamplingStrategy()) {
    case 1:
      local_problem_op.assemble_matrix(locprob_system_matrix);
      local_problem_op.assemble_local_RHS(e, local_problem_rhs);
      break;
    case 2:
      if (coarse_index < 0)
        DUNE_THROW(Dune::InvalidStateException, "Invalid coarse index: coarse_index < 0");
      local_problem_op.assemble_matrix(locprob_system_matrix, subgrid_list_.getCoarseNodeVector(coarse_index));
      local_problem_op.assemble_local_RHS(e, subgrid_list_.getCoarseNodeVector(coarse_index),
                                          specifier_.getOversamplingStrategy(), local_problem_rhs);
      break;
    case 3:
      if (coarse_index < 0)
        DUNE_THROW(Dune::InvalidStateException, "Invalid coarse index: coarse_index < 0");
      if (DSC_CONFIG_GET("rigorous_msfem.oversampling_strategy", "Clement") == "Clement") {
        local_problem_op.assemble_matrix(locprob_system_matrix);
      } else {
        local_problem_op.assemble_matrix(locprob_system_matrix, subgrid_list_.getCoarseNodeVector(coarse_index));
      }
      local_problem_op.assemble_local_RHS(e, subgrid_list_.getCoarseNodeVector(coarse_index),
                                          specifier_.getOversamplingStrategy(), local_problem_rhs);
      break;
    default:
      DUNE_THROW(Dune::Fem::ParameterInvalid, "Oversampling Strategy must be 1, 2 or 3.");
  }

  //! boundary treatment:
  for (const auto& subgridEntity : subDiscreteFunctionSpace) {
    DSFe::LocalMatrixProxy<LocProbLinearOperatorTypeType> localMatrix(locprob_system_matrix, subgridEntity,
                                                                      subgridEntity);

    const auto& lagrangePointSet = subDiscreteFunctionSpace.lagrangePointSet(subgridEntity);
    auto rhsLocal = local_problem_rhs.localFunction(subgridEntity);

    for (const auto& subgridIntersection : DSC::intersectionRange(subDiscreteFunctionSpace.gridPart(), subgridEntity)) {
      // if there is a neighbor entity
      if (subgridIntersection.boundary()) {
        const int face = subgridIntersection.indexInInside();
        for (const auto& lp : DSC::lagrangePointSetRange(lagrangePointSet, face)) {
          // zero boundary condition for 'cell problems':
          // set unit row in matrix for any boundary dof ...
          localMatrix.unitRow(lp);
          // ... and set respective rhs dof to zero
          rhsLocal[lp] = 0;
        }
      }
    }
  }

  if (!(local_problem_rhs.dofsValid())) {
    DUNE_THROW(Dune::InvalidStateException, "Local MsFEM Problem RHS invalid.");
  }

  // NOTE:
  // is the right hand side of the local MsFEM problem equal to zero or almost identical to zero?
  // if yes, the solution of the local MsFEM problem is also identical to zero. The solver is getting a problem with
  // this situation, which is why we do not solve local msfem problems for zero-right-hand-side, since we already know
  // the result.
  if (local_problem_op.normRHS(local_problem_rhs) < /*1e-06*/ 1e-30) {
    local_problem_solution.clear();
    DSC_LOG_ERROR << "Local MsFEM problem with solution zero." << std::endl;
    return;
  }

  const auto locprob_fem_biCGStab = make_inverse_operator(locprob_system_matrix);

  const bool clement = specifier_.getOversamplingStrategy() == 3
                           ? DSC_CONFIG_GET("rigorous_msfem.oversampling_strategy", "Clement") == "Clement"
                           : false;

  if (clement) {
    HostDiscreteFunctionType zero("zero", specifier_.coarseSpace());
    zero.clear();
    const double dummy = 12345.67890;
    double solverEps = 1e-2;
    int maxIterations = 1000;

    // we want to solve the local problem with the constraint that the weighted Clement interpoltion
    // of the local problem solution is zero

    // implementation of a weighted Clement interpolation operator for our purpose:
    WeightedClementOperatorType clement_interpolation_op(subDiscreteFunctionSpace, specifier_.coarseSpace(),
                                                         subgrid_list_.getCoarseNodeVector(coarse_index),
                                                         *coarse_basis_, *global_id_to_internal_id_, specifier_);
    //! NOTE TODO: implementation is not yet optimal, because the weighted Clement maps a function
    //! defined on the local subgrid to a function defined on the whole(!) coarse space.
    //! It would be better to implement a mapping to a localized coarse space, since
    //! the uzawa solver must treat ALL coarse grid nodes (expensive and worse convergence).

    // saddle point problem solver:
    typedef UzawaInverseOp<SubDiscreteFunctionType, HostDiscreteFunctionType, InverseLocProbLinearOperatorTypeType,
                           WeightedClementOperatorType> InverseUzawaOperatorType;

    HostDiscreteFunctionType lagrange_multiplier("lagrange multiplier", specifier_.coarseSpace());
    lagrange_multiplier.clear();

    // create inverse operator
    // saddle point problem solver with uzawa algorithm:
    {
      DSC::Profiler::ScopedTiming st("uzawa");
      InverseUzawaOperatorType uzawa(*locprob_fem_biCGStab, clement_interpolation_op, dummy, solverEps, maxIterations,
                                     true);
      uzawa(local_problem_rhs, zero /*interpolation is zero*/, local_problem_solution, lagrange_multiplier);
    }
  } else {
    locprob_fem_biCGStab->apply(local_problem_rhs, local_problem_solution);
  }

  if (!(local_problem_solution.dofsValid())) {
    DUNE_THROW(Dune::InvalidStateException, "Current solution of the local msfem problem invalid!");
  }
} // solvelocalproblem

void MsFEMLocalProblemSolver::subgrid_to_hostrid_function(const SubDiscreteFunctionType& sub_func,
                                                          HostDiscreteFunctionType& host_func) const {
  host_func.clear();

  const auto& subDiscreteFunctionSpace = sub_func.space();
  const auto& subGrid = subDiscreteFunctionSpace.grid();
  for (const auto& sub_entity : subDiscreteFunctionSpace) {
    auto host_entity_pointer = subGrid.getHostEntity<0>(sub_entity);
    const auto& host_entity = *host_entity_pointer;

    auto sub_loc_value = sub_func.localFunction(sub_entity);
    auto host_loc_value = host_func.localFunction(host_entity);

    const auto numBaseFunctions = sub_loc_value.basisFunctionSet().size();
    for (unsigned int i = 0; i < numBaseFunctions; ++i) {
      host_loc_value[i] = sub_loc_value[i];
    }
  }
} // subgrid_to_hostrid_function

void MsFEMLocalProblemSolver::output_local_solution(const int coarse_index, const int which,
                                                    const HostDiscreteFunctionType& host_local_solution) const {
  if (!DSC_CONFIG_GET("lod.local_solution_vtk_output", true))
    return;
  typedef tuple<const HostDiscreteFunctionType*> IOTupleType;
  typedef Fem::DataOutput<HostGridType, IOTupleType> DataOutputType;

  // general output parameters
  LocalProblemDataOutputParameters outputparam;
  // --------- data output local solution --------------

  // create and initialize output class
  IOTupleType local_solution_series(&host_local_solution);

  const std::string ls_name_s = (boost::format("local_problem_solution_e%d_%d") % which % coarse_index).str();
  // std::cout << "Write local corrector in vtk format to: " << ls_name_s << std::endl;
  outputparam.set_prefix(ls_name_s);
  DataOutputType localsol_dataoutput(hostDiscreteFunctionSpace_.gridPart().grid(), local_solution_series, outputparam);
  localsol_dataoutput.writeData(1.0 /*dummy*/, (boost::format("local-problem-solution-%d") % which).str());
}

void MsFEMLocalProblemSolver::output_local_solution(const int coarseIndex, const int which,
                                                    const SubDiscreteFunctionType& solution) const {
  HostDiscreteFunctionType hostSolution("Local Solution", hostDiscreteFunctionSpace_);
  hostSolution.clear();
  subgrid_to_hostrid_function(solution, hostSolution);
  output_local_solution(coarseIndex, which, hostSolution);
}

void MsFEMLocalProblemSolver::assembleAndSolveAll(bool /*silent*/) {
  const bool uzawa = DSC_CONFIG_GET("rigorous_msfem.uzawa_solver", false);
  const bool clement = (DSC_CONFIG_GET("rigorous_msfem.oversampling_strategy", "Clement") == "Clement");
  if ((!uzawa) && (specifier_.getOversamplingStrategy() == 3) && clement) {
    DUNE_THROW(NotImplemented, "caller needs to instantiate LodLocalProblemSolver and call assemble_all_clement_lod on it");
  }
  if (uzawa && !(specifier_.simplexCoarseGrid())) {
    DUNE_THROW(NotImplemented, "Uzawa-solver and non-simplex grid have not been tested together, yet!");
  }

  // number of coarse grid entities (of codim 0).
  const auto coarseGridSize = specifier_.getNumOfCoarseEntities();

  DSC_LOG_INFO << "in method 'assemble_all': coarseGridSize = " << coarseGridSize << std::endl;
  DSC_PROFILER.startTiming("msfem.localproblemsolver.assemble_all");

  threadIterators_.update();
  // we want to determine minimum, average and maxiumum time for solving a local msfem problem in the current method
  DSC::MinMaxAvg<double> cell_time;

  // stupid way to pre-init the discretefunctions spaces in the df_io backend inside a serial section
  const auto& coarseSpace = specifier_.coarseSpace();
  for (const auto& coarseEntity : coarseSpace) {
    LocalSolutionManager(coarseEntity, subgrid_list_, specifier_).solutionsWereLoaded();
  }

  #ifdef _OPENMP
  #pragma omp parallel
  #endif
  {
  const auto& coarseGridLeafIndexSet = coarseSpace.gridPart().grid().leafIndexSet();
  for (const auto& coarseEntity : threadIterators_) {
    const int coarse_index = coarseGridLeafIndexSet.index(coarseEntity);
    
    DSC_LOG_INFO << "-------------------------" << std::endl << "Coarse index " << coarse_index << std::endl;
    DSC_PROFILER.startTiming("none.local_problem_solution");
    LocalSolutionManager localSolutionManager(coarseEntity, subgrid_list_, specifier_);
    // solve the problems
    solveAllLocalProblems(coarseEntity, localSolutionManager.getLocalSolutions());

    cell_time(DSC_PROFILER.stopTiming("none.local_problem_solution") / 1000.f);
    DSC_PROFILER.resetTiming("none.local_problem_solution");
    localSolutionManager.saveLocalSolutions();
  } // for
  } // omp region

  //! @todo The following debug-output is wrong (number of local problems may be different)
  const auto total_time = DSC_PROFILER.stopTiming("msfem.localproblemsolver.assemble_all") / 1000.f;
  DSC_LOG_INFO << std::endl;
  DSC_LOG_INFO << "In method: assemble_all." << std::endl << std::endl;
  DSC_LOG_INFO << "MsFEM problems solved for " << coarseGridSize << " coarse grid entities." << std::endl;
  DSC_LOG_INFO << CommonTraits::GridType::dimension * coarseGridSize << " local MsFEM problems solved in total." << std::endl;
  DSC_LOG_INFO << "Minimum time for solving a local problem = " << cell_time.min() << "s." << std::endl;
  DSC_LOG_INFO << "Maximum time for solving a localproblem = " << cell_time.max() << "s." << std::endl;
  DSC_LOG_INFO << "Average time for solving a localproblem = " << cell_time.average() << "s." << std::endl;
  DSC_LOG_INFO << "Total time for computing and saving the localproblems = " << total_time << "s," << std::endl
               << std::endl;
} // assembleAndSolveAll

} // namespace MsFEM {
} // namespace Multiscale {
} // namespace Dune {
