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

std::unique_ptr<LocalProblemSolver::InverseLocProbLinearOperatorTypeType>
LocalProblemSolver::make_inverse_operator(LocalProblemSolver::LocProbLinearOperatorTypeType& problem_matrix) {
  const auto solver =
      Dune::Multiscale::Problem::getModelData()->symmetricDiffusion() ? std::string("cg") : std::string("bcgs");
  return DSC::make_unique<InverseLocProbLinearOperatorTypeType>(
      problem_matrix, 1e-8, 1e-8, 20000, DSC_CONFIG_GET("localproblemsolver_verbose", false), solver,
      DSC_CONFIG_GET("preconditioner_type", std::string("sor")), 1);
}

LocalProblemSolver::LocalProblemSolver(const MacroMicroGridSpecifierType& specifier,
                                                 LocalGridList& subgrid_list,
                                                 const DiffusionOperatorType& diffusion_operator)
  : diffusion_(diffusion_operator)
  , specifier_(specifier)
  , subgrid_list_(subgrid_list)
  #ifdef ENABLE_LOD_ONLY_CODE
  , ids_relevant_basis_functions_for_subgrid_(nullptr)
  , inverse_of_L1_norm_coarse_basis_funcs_(nullptr)
  , coarse_basis_(nullptr)
  , global_id_to_internal_id_(nullptr)
  , neumann_bc_(nullptr)
  , dirichlet_extension_(nullptr)
  #endif // ENABLE_LOD_ONLY_CODE
{}

LocalProblemSolver::LocalProblemSolver(const MacroMicroGridSpecifierType& specifier,
    LocalGridList& subgrid_list, std::vector<std::vector<std::size_t>>& /*ids_basis_functions_in_subgrid*/,
    std::vector<double>& /*inverse_of_L1_norm_coarse_basis_funcs*/, const DiffusionOperatorType& diffusion_operator,
    const CoarseBasisFunctionListType& /*coarse_basis*/, const std::map<std::size_t, std::size_t>& /*global_id_to_internal_id*/,
    const NeumannBoundaryType& /*neumann_bc*/, const LocalGridDiscreteFunctionType& /*dirichlet_extension*/)
  : diffusion_(diffusion_operator)
  , specifier_(specifier)
  , subgrid_list_(subgrid_list)
  #ifdef ENABLE_LOD_ONLY_CODE
  , ids_relevant_basis_functions_for_subgrid_(&ids_basis_functions_in_subgrid)
  ,
  // ids of the coarse grid basis functions in the interior of the subgrid
  inverse_of_L1_norm_coarse_basis_funcs_(&inverse_of_L1_norm_coarse_basis_funcs)
  , coarse_basis_(&coarse_basis)
  , global_id_to_internal_id_(&global_id_to_internal_id)
  , neumann_bc_(&neumann_bc)
  , dirichlet_extension_(&dirichlet_extension)
  #endif // ENABLE_LOD_ONLY_CODE
{}

/** Solve all local MsFEM problems for one coarse entity at once.
*
*
*/
void LocalProblemSolver::solveAllLocalProblems(const CoarseEntityType& coarseCell,
                                                    LocalGridDiscreteFunctionVectorType& allLocalSolutions) const {
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
  LocalGridDiscreteFunctionVectorType allLocalRHS(allLocalSolutions.size());
  for (auto& it : allLocalRHS)
    it = DSC::make_unique<LocalGridDiscreteFunctionType>("rhs of local MsFEM problem", subDiscreteFunctionSpace);

  if (DSC_CONFIG_GET("msfem.oversampling_strategy", 1) != 1)
    DUNE_THROW(Fem::ParameterInvalid, "Oversampling Strategy must be 1 at the moment");

  localProblemOperator.assemble_matrix(locProbSysMatrix);
  localProblemOperator.assembleAllLocalRHS(coarseCell, specifier_, allLocalRHS);

  // set dirichlet dofs to zero
  Stuff::GridboundaryAllDirichlet<LocalGridType::LeafGridView::Intersection> boundaryInfo;
  DirichletConstraints<LocalGridDiscreteFunctionSpaceType> constraints(boundaryInfo, subDiscreteFunctionSpace);
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

void LocalProblemSolver::output_local_solution(const int coarse_index, const int which,
                                                    const LocalGridDiscreteFunctionType& host_local_solution) const {
  if (!DSC_CONFIG_GET("lod.local_solution_vtk_output", true))
    return;
  typedef tuple<const LocalGridDiscreteFunctionType*> IOTupleType;
  typedef Fem::DataOutput<LocalGridType, IOTupleType> DataOutputType;

  // general output parameters
  LocalProblemDataOutputParameters outputparam;
  // --------- data output local solution --------------

  // create and initialize output class
  IOTupleType local_solution_series(&host_local_solution);

  const std::string ls_name_s = (boost::format("local_problem_solution_e%d_%d") % which % coarse_index).str();
  // std::cout << "Write local corrector in vtk format to: " << ls_name_s << std::endl;
  outputparam.set_prefix(ls_name_s);
  DataOutputType localsol_dataoutput(host_local_solution.space().grid(), local_solution_series, outputparam);
  localsol_dataoutput.writeData(1.0 /*dummy*/, (boost::format("local-problem-solution-%d") % which).str());
}

void LocalProblemSolver::assembleAndSolveAll(bool /*verbose*/) {
  static const int dimension = CommonTraits::GridType::dimension;

  JacobianRangeType unitVectors[dimension];
  for (int i = 0; i < dimension; ++i)
    for (int j = 0; j < dimension; ++j) {
      if (i == j) {
        unitVectors[i][0][j] = 1.0;
      } else {
        unitVectors[i][0][j] = 0.0;
      }
    }

  // number of coarse grid entities (of codim 0).
  const auto coarseGridSize = specifier_.coarseSpace().grid().size(0);
  if (Dune::Fem::MPIManager::size()>0)
    DSC_LOG_INFO << "Rank " << Dune::Fem::MPIManager::rank()
                 << " will solve local problems for " << coarseGridSize
                 << " coarse entities!" << std::endl;
    else {
    DSC_LOG_INFO << "Will solve local problems for " << coarseGridSize
                 << " coarse entities!" << std::endl;
  }
  DSC_PROFILER.startTiming("msfem.localproblemsolver.assemble_all");

  // we want to determine minimum, average and maxiumum time for solving a local msfem problem in the current method
  DSC::MinMaxAvg<double> cell_time;
  DSC::MinMaxAvg<double> saveTime;

  const auto& coarseSpace = specifier_.coarseSpace();
  const auto& coarseGridLeafIndexSet = coarseSpace.gridPart().grid().leafIndexSet();
  for (const auto& coarseEntity : coarseSpace) {
    const int coarse_index = coarseGridLeafIndexSet.index(coarseEntity);

    DSC_LOG_INFO << "-------------------------" << std::endl << "Coarse index " << coarse_index << std::endl;

    DSC_PROFILER.startTiming("none.saveLocalProblemsOnCell");

    const bool uzawa = DSC_CONFIG_GET("rigorous_msfem.uzawa_solver", false);
    const bool clement = (DSC_CONFIG_GET("rigorous_msfem.oversampling_strategy", "Clement") == "Clement");
    if ((!uzawa) && (DSC_CONFIG_GET("msfem.oversampling_strategy", 1) == 3) && clement) {
      DUNE_THROW(InvalidStateException, "broken lod only code");
    #ifdef ENABLE_LOD_ONLY_CODE
      //! only for dimension 2 and simplex grid!

      // preprocessing step to assemble the two relevant system matrices:
      // one for the corrector problem without contraints
      // and the second of the low dimensional lagrange multiplier
      // (describing the inverse of the schur complement)
      // the preprocessing step takes the main costs for solving the corrector problems,
      // however, it is the same for all correctors on the same subgrid
      // -----------------------------------------------------------------------------------------------------
      //! the matrix in our linear system of equations
      // in the non-linear case, it is the matrix for each iteration step
      LocProbLinearOperatorTypeType locprob_system_matrix("Local Problem System Matrix", subDiscreteFunctionSpace,
                                                          subDiscreteFunctionSpace);

      const auto number_of_relevant_coarse_nodes_for_subgrid =
          (*ids_relevant_basis_functions_for_subgrid_)[coarse_index].size();
      MatrixType lagrange_multiplier_system_matrix(number_of_relevant_coarse_nodes_for_subgrid,
                                                   number_of_relevant_coarse_nodes_for_subgrid);

      std::cout << "Start preprocessing for subgrid " << coarse_index << "." << std::endl;
      std::cout << "There are " << number_of_relevant_coarse_nodes_for_subgrid << " preprocessing problems to solve."
                << std::endl;
      preprocess_corrector_problems(coarse_index, locprob_system_matrix, lagrange_multiplier_system_matrix);
      std::cout << "Preprocessing done." << std::endl;
      // -----------------------------------------------------------------------------------------------------

      LocalGridDiscreteFunctionType local_problem_solution_0(name_local_solution, subDiscreteFunctionSpace);
      local_problem_solution_0.clear();

      LocalGridDiscreteFunctionType local_problem_solution_1(name_local_solution, subDiscreteFunctionSpace);
      local_problem_solution_1.clear();

      // 'solve' methods requires the pre-processing step (that is the same for both directions e_0 and e_1)
      // (this at least halfes the computational complexity)
      solve_corrector_problem_lod(unitVectors[0], locprob_system_matrix, lagrange_multiplier_system_matrix,
                                  local_problem_solution_0, coarse_index);
      solve_corrector_problem_lod(unitVectors[1], locprob_system_matrix, lagrange_multiplier_system_matrix,
                                  local_problem_solution_1, coarse_index);

      assert(local_problem_solution_0.dofsValid());
      assert(local_problem_solution_1.dofsValid());

      const std::string locprob_solution_location =
          (boost::format("local_problems/_localProblemSolutions_%d") % coarseId).str();
      DiscreteFunctionWriter dfw(locprob_solution_location);
      dfw.append(local_problem_solution_0);
      dfw.append(local_problem_solution_1);

      std::vector<std::size_t> indices;
      coarseSpace.mapper().map(coarseEntity, indices);

      // was the dirichlet (resp. neumann) boundary corrector for this element already assembled?
      bool dirichlet_boundary_corrector_assembled = false;
      bool neumann_boundary_corrector_assembled = false;
      // assemble Dirichlet and Neumann boundary correctors
      for (const auto& intersection : DSC::intersectionRange(coarseSpace.gridPart(), coarseEntity)) {
        if ((!dirichlet_extension_) || (!neumann_bc_))
          continue;

        bool solve_for_dirichlet_corrector = false;
        if (!dirichlet_boundary_corrector_assembled) {
          if (intersection.boundary() && (intersection.boundaryId() == 1)) {
            solve_for_dirichlet_corrector = true;
          } else {
            const auto& lagrangePointSet = coarseSpace.lagrangePointSet(coarseEntity);

            const int face = intersection.indexInInside();
            for (const auto& lp : DSC::lagrangePointSetRange(lagrangePointSet, face)) {
              if (specifier_.is_coarse_dirichlet_node(indices[lp])) {
                solve_for_dirichlet_corrector = true;
              }
            }
            if ((!solve_for_dirichlet_corrector) && (!intersection.boundary())) {
              continue;
            }
          }
        } else {
          solve_for_dirichlet_corrector = false;
        }

        // Dirichlet boundary corrector:
        if (solve_for_dirichlet_corrector) {
          const std::string name_dirichlet_corrector =
              (boost::format("Dirichlet Boundary Corrector %d") % coarseId).str();
          LocalGridDiscreteFunctionType dirichlet_boundary_corrector(name_dirichlet_corrector, subDiscreteFunctionSpace);
          dirichlet_boundary_corrector.clear();

          // also requires the pre-processing step:
          std::cout << "Solve Dirichlet boundary corrector problem for subgrid " << coarse_index << std::endl;
          solve_dirichlet_corrector_problem_lod(locprob_system_matrix, lagrange_multiplier_system_matrix,
                                                dirichlet_boundary_corrector, coarse_index);

          assert(dirichlet_boundary_corrector.dofsValid());

          const std::string dirichlet_corrector_location =
              (boost::format("local_problems/_dirichletBoundaryCorrector_%d") % coarseId).str();
          DiscreteFunctionWriter dfw_dirichlet(dirichlet_corrector_location);
          dfw_dirichlet.append(dirichlet_boundary_corrector);
          dirichlet_boundary_corrector_assembled = true;

          if (DSC_CONFIG_GET("lod.localproblem_vtkoutput", true)) {
            output_local_solution(coarse_index, 2, dirichlet_boundary_corrector); // 2 stands for Dirichlet
          }
        }

        // Neumann boundary corrector:
        if (intersection.boundary() && (intersection.boundaryId() == 2) && (!neumann_boundary_corrector_assembled)) {
          const std::string name_neumann_corrector = (boost::format("Neumann Boundary Corrector %d") % coarseId).str();
          LocalGridDiscreteFunctionType neumann_boundary_corrector(name_neumann_corrector, subDiscreteFunctionSpace);
          neumann_boundary_corrector.clear();
          std::cout << "Solve Neumann boundary corrector problem for subgrid " << coarse_index << std::endl;

          // also requires the pre-processing step:
          solve_neumann_corrector_problem_lod(locprob_system_matrix, lagrange_multiplier_system_matrix,
                                              neumann_boundary_corrector, coarse_index);

          assert(neumann_boundary_corrector.dofsValid());

          const std::string neumann_corrector_location =
              (boost::format("local_problems/_neumannBoundaryCorrector_%d") % coarseId).str();
          DiscreteFunctionWriter dfw_neumann(neumann_corrector_location);
          dfw_neumann.append(neumann_boundary_corrector);
          neumann_boundary_corrector_assembled = true;

          if (DSC_CONFIG_GET("lod.localproblem_vtkoutput", true)) {
            output_local_solution(coarse_index, 3, neumann_boundary_corrector); // 3 stands for Neumann
          }
        }
      }

      if (DSC_CONFIG_GET("lod.localproblem_vtkoutput", true)) {
        output_local_solution(coarse_index, 0, local_problem_solution_0);
        output_local_solution(coarse_index, 1, local_problem_solution_1);
      }
    #endif // ENABLE_LOD_ONLY_CODE
    } else if (uzawa && !(specifier_.simplexCoarseGrid())) {
      DUNE_THROW(NotImplemented, "Uzawa-solver and non-simplex grid have not been tested together, yet!");
    } else {
      // take time
      DSC_PROFILER.startTiming("none.local_problem_solution");
      LocalSolutionManager localSolutionManager(coarseEntity, subgrid_list_, specifier_);

      // solve the problems
      solveAllLocalProblems(coarseEntity, localSolutionManager.getLocalSolutions());
      // min/max time
      cell_time(DSC_PROFILER.stopTiming("none.local_problem_solution") / 1000.f);
      DSC_PROFILER.resetTiming("none.local_problem_solution");

      // save the local solutions to disk
      DSC_PROFILER.startTiming("none.saveLocalProblemSolution");
      localSolutionManager.save();
      saveTime(DSC_PROFILER.stopTiming("none.saveLocalProblemSolution") / 1000.f);
      DSC_PROFILER.resetTiming("none.saveLocalProblemSolution");
    }

    DSC_LOG_INFO << "Total time for solving and saving all local problems for the current subgrid: "
                 << DSC_PROFILER.stopTiming("none.saveLocalProblemsOnCell") / 1000.f << "s"
                 << std::endl << std::endl;
    DSC_PROFILER.resetTiming("none.saveLocalProblemsOnCell");
  } // for

  //! @todo The following debug-output is wrong (number of local problems may be different)
  const auto total_time = DSC_PROFILER.stopTiming("msfem.localproblemsolver.assemble_all") / 1000.f;
  DSC_LOG_INFO << std::endl;
  DSC_LOG_INFO << "Local problems solved for " << coarseGridSize << " coarse grid entities." << std::endl;
  DSC_LOG_INFO << "Minimum time for solving a local problem = " << cell_time.min() << "s." << std::endl;
  DSC_LOG_INFO << "Maximum time for solving a localproblem = " << cell_time.max() << "s." << std::endl;
  DSC_LOG_INFO << "Average time for solving a localproblem = " << cell_time.average() << "s." << std::endl;
  DSC_LOG_INFO << "Minimum time for saving a local problem = " << saveTime.min() << "s." << std::endl;
  DSC_LOG_INFO << "Maximum time for saving a localproblem = " << saveTime.max() << "s." << std::endl;
  DSC_LOG_INFO << "Average time for saving a localproblem = " << saveTime.average() << "s." << std::endl;
  DSC_LOG_INFO << "Total time for computing and saving the localproblems = " << total_time << "s," << std::endl
               << std::endl;
} // assemble_all

#ifdef ENABLE_LOD_ONLY_CODE
void LocalProblemSolver::solvelocalproblem(JacobianRangeType& e, LocalGridDiscreteFunctionType& local_problem_solution,
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
  LocalGridDiscreteFunctionType local_problem_rhs("rhs of local MsFEM problem", subDiscreteFunctionSpace);
  local_problem_rhs.clear();
  if (coarse_index < 0)
    DUNE_THROW(Dune::InvalidStateException, "Invalid coarse index: coarse_index < 0");
  const auto oversampling_strategy = DSC_CONFIG_GET("msfem.oversampling_strategy", 1);
  switch (oversampling_strategy) {
    case 1:
      local_problem_op.assemble_matrix(locprob_system_matrix);
      local_problem_op.assemble_local_RHS(e, local_problem_rhs);
      break;
    case 2:
      local_problem_op.assemble_matrix(locprob_system_matrix, subgrid_list_.getCoarseNodeVector(coarse_index));
      local_problem_op.assemble_local_RHS(e, subgrid_list_.getCoarseNodeVector(coarse_index),
                                          oversampling_strategy, local_problem_rhs);
      break;
    case 3:
      if (DSC_CONFIG_GET("rigorous_msfem.oversampling_strategy", "Clement") == "Clement") {
        local_problem_op.assemble_matrix(locprob_system_matrix);
      } else {
        local_problem_op.assemble_matrix(locprob_system_matrix, subgrid_list_.getCoarseNodeVector(coarse_index));
      }
      local_problem_op.assemble_local_RHS(e, subgrid_list_.getCoarseNodeVector(coarse_index),
                                          oversampling_strategy, local_problem_rhs);
      break;
    default:
      DUNE_THROW(Dune::Fem::ParameterInvalid, "Oversampling Strategy must be 1, 2 or 3.");
  }

  //! boundary treatment:
  for (const auto& local_entity : subDiscreteFunctionSpace) {
    DSFe::LocalMatrixProxy<LocProbLinearOperatorTypeType> localMatrix(locprob_system_matrix, local_entity,
                                                                      local_entity);

    const auto& lagrangePointSet = subDiscreteFunctionSpace.lagrangePointSet(local_entity);
    auto rhsLocal = local_problem_rhs.localFunction(local_entity);

    for (const auto& subgridIntersection : DSC::intersectionRange(subDiscreteFunctionSpace.gridPart(), local_entity)) {
      if (subgridIntersection.boundary()) {
        const int face = subgridIntersection.indexInInside();
        for (const auto& lp : DSC::lagrangePointSetRange(lagrangePointSet, face)) {
          localMatrix.unitRow(lp);
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

  const bool clement = DSC_CONFIG_GET("msfem.oversampling_strategy", 1) == 3
                           ? DSC_CONFIG_GET("rigorous_msfem.oversampling_strategy", "Clement") == "Clement"
                           : false;

  if (clement) {
    CommonTraits::DiscreteFunctionType zero("zero", specifier_.coarseSpace());
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
    typedef UzawaInverseOp<LocalGridDiscreteFunctionType, LocalGridDiscreteFunctionType, InverseLocProbLinearOperatorTypeType,
                           WeightedClementOperatorType> InverseUzawaOperatorType;

    CommonTraits::DiscreteFunctionType lagrange_multiplier("lagrange multiplier", specifier_.coarseSpace());
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

// assemble the two relevant system matrices: one for the corrector problem without contraints
// and the second of the low dimensional lagrange multiplier (describing the inverse of the schur complement)
void LocalProblemSolver::preprocess_corrector_problems(const int coarse_index,
                                                            LocProbLinearOperatorTypeType& locprob_system_matrix,
                                                            MatrixType& lm_system_matrix) const {
  auto subGridPart = subgrid_list_.gridPart(coarse_index);
  const LocalGridDiscreteFunctionSpaceType subDiscreteFunctionSpace(subGridPart);

  //! define the discrete (elliptic) local MsFEM problem operator
  // ( effect of the discretized differential operator on a certain discrete function )
  LocalProblemOperator local_problem_op(subDiscreteFunctionSpace, diffusion_);

  const auto& subGrid = subDiscreteFunctionSpace.grid();

  // NOTE:
  // is the right hand side of the local MsFEM problem equal to zero or almost identical to zero?
  // if yes, the solution of the local MsFEM problem is also identical to zero. The solver is getting a problem with
  // this situation, which is why we do not solve local msfem problems for zero-right-hand-side, since we already know
  // the result.

  switch (DSC_CONFIG_GET("msfem.oversampling_strategy", 1)) {
    case 3:
      break;
    default:
      DUNE_THROW(Dune::InvalidStateException,
                 "method 'preprocess_corrector_problems' can be only used in combination with the LOD.");
  }

  bool clement = (DSC_CONFIG_GET("rigorous_msfem.oversampling_strategy", "Clement") == "Clement");

  if (!clement)
    DUNE_THROW(Dune::InvalidStateException, "method 'preprocess_corrector_problems' can be only used in combination "
                                            "with the LOD and Clement interpolation.");
  local_problem_op.assemble_matrix(locprob_system_matrix);

  //! boundary treatment:
  const auto& localGridPart = hostDiscreteFunctionSpace_.gridPart();
  for (const auto& localgrid_entity : subDiscreteFunctionSpace) {
//    auto host_entity_pointer = subGrid.getLocalEntity<0>(localgrid_entity);
    const auto& host_entity = localgrid_entity;//*host_entity_pointer;

    DSFe::LocalMatrixProxy<LocProbLinearOperatorTypeType> local_matrix(locprob_system_matrix, localgrid_entity,
                                                                       localgrid_entity);

    for (const auto& intersection : DSC::intersectionRange(localGridPart, host_entity)) {
      if (intersection.neighbor()) {
        // check if the neighbor entity is in the subgrid
        const auto neighborLocalEntityPointer = intersection.outside();
        const auto& neighborLocalEntity = *neighborLocalEntityPointer;
        if (subGrid.contains<0>(neighborLocalEntity)) {
          continue;
        }
      } else if (intersection.boundaryId() != 1) // check if the intersection is part of the Dirichlet boundary
      {
        continue;
      }

      const auto face = intersection.indexInInside();
      const auto& lagrangePointSet = subDiscreteFunctionSpace.lagrangePointSet(localgrid_entity);
      for (const auto& lp : DSC::lagrangePointSetRange(lagrangePointSet, face))
        local_matrix.unitRow(lp);
    }
  }
  // ----------------------------------------------------------------------------------------------------

  //! Essential pre-processing step
  // For each coarse node j in the subgrid (local internal numbering, i.e. 0 <= j < M_subgrid), solve
  // for b_h_j with S_h b_h_j = C_h^T e_j, where C_h describes the algebraic version of the weighted
  // Clement interpolation operator
  // ----------------------------------------------------------------------------------------------------
  auto number_of_relevant_coarse_nodes_for_subgrid = (*ids_relevant_basis_functions_for_subgrid_)[coarse_index].size();

  assert(number_of_relevant_coarse_nodes_for_subgrid);
  std::vector<std::unique_ptr<LocalGridDiscreteFunctionType>> b_h(number_of_relevant_coarse_nodes_for_subgrid);
  std::vector<std::unique_ptr<LocalGridDiscreteFunctionType>> rhs_Chj(number_of_relevant_coarse_nodes_for_subgrid);
  for (auto j : DSC::valueRange(number_of_relevant_coarse_nodes_for_subgrid)) {
    b_h[j] = DSC::make_unique<LocalGridDiscreteFunctionType>("q_h", subDiscreteFunctionSpace);
    rhs_Chj[j] = DSC::make_unique<LocalGridDiscreteFunctionType>("rhs_Chj_h", subDiscreteFunctionSpace);
    b_h[j]->clear();
    rhs_Chj[j]->clear();
  }

  // clement_weight_j ( \psi_i, \Psi_j ), where
  // clement_weight_j = (*inverse_of_L1_norm_coarse_basis_funcs_)[interior_basis_func_id]
  local_problem_op.assemble_local_RHS_lg_problems_all((*coarse_basis_), (*inverse_of_L1_norm_coarse_basis_funcs_),
                                                      (*ids_relevant_basis_functions_for_subgrid_)[coarse_index],
                                                      rhs_Chj);
  // get the global id of all interior coarse basis functions (subgrid id, local id) -> (global interior id)
  // zero boundary condition for 'rhs_Chj[j]' (except on the Neumann boundary part):
  // set Dirichlet Boundary to zero
  for (const auto& localgrid_entity : subDiscreteFunctionSpace) {
    auto host_entity_pointer = subGrid.getLocalEntity<0>(localgrid_entity);
    const auto& host_entity = *host_entity_pointer;

    HostIntersectionIterator iit = localGridPart.ibegin(host_entity);
    const HostIntersectionIterator endiit = localGridPart.iend(host_entity);
    for (; iit != endiit; ++iit) {
      if (iit->neighbor()) // if there is a neighbor entity
      {
        // check if the neighbor entity is in the subgrid
        const auto neighborLocalEntityPointer = iit->outside();
        const auto& neighborLocalEntity = *neighborLocalEntityPointer;

        if (subGrid.contains<0>(neighborLocalEntity)) {
          continue;
        }
      } else if (iit->boundaryId() != 1) // check if the intersection is part of the Dirichlet boundary
      {
        continue;
      }

      // SubLocalFunctionType rhsLocal = rhs_Chj[j]->localFunction(localgrid_entity);
      const auto& lagrangePointSet = subDiscreteFunctionSpace.lagrangePointSet(localgrid_entity);
      const auto face = (*iit).indexInInside();
      for (const auto& lp : DSC::lagrangePointSetRange(lagrangePointSet, face))
        for (auto j : DSC::valueRange(number_of_relevant_coarse_nodes_for_subgrid))
          ((rhs_Chj[j])->localFunction(localgrid_entity))[lp] = 0;
    }
  }

  const InverseLocProbLinearOperatorTypeType locprob_inverse_system_matrix(
      locprob_system_matrix, 1e-8, 1e-8, 20000, DSC_CONFIG_GET("lod.local_problem_solver_verbose", false),
      DSC_CONFIG_GET("lod.local_solver", "bcgs"), DSC_CONFIG_GET("preconditioner_type", std::string("sor")));
  // solve the pre-processing problems:
  for (auto j : DSC::valueRange(number_of_relevant_coarse_nodes_for_subgrid))
    locprob_inverse_system_matrix.apply(*(rhs_Chj[j]), *(b_h[j]));

  // ----------------------------------------------------------------------------------------------------

  //! Assemble the (low dimensional) system matrix for the final lagrange multiplier problem
  // (matrix of size 'number_of_relevant_coarse_nodes_for_subgrid')
  // ----------------------------------------------------------------------------------------------------

  // (stiffness) matrix for the lagrange multiplier (lm) problem
  for (size_t i = 0; i != lm_system_matrix.N(); ++i)   // rows
    for (size_t j = 0; j != lm_system_matrix.M(); ++j) // colums
      lm_system_matrix[i][j] = 0.0;

  // matrix with entries M[i][j] = weight_i ( b_h[j], coarse_basis_func[i] )_L2(\Omega)
  // 'i = row' and 'j = column'

  for (const auto& localgrid_entity : subDiscreteFunctionSpace) {
    const auto& sg_geometry = localgrid_entity.geometry();

    auto host_entity_pointer = subGrid.getLocalEntity<0>(localgrid_entity);
    const auto& host_entity = *host_entity_pointer;

    // exact for polynomials of degree 2:
    const auto sg_quadrature = DSFe::make_quadrature(localgrid_entity, subDiscreteFunctionSpace);
    const auto quadrature = DSFe::make_quadrature(host_entity, hostDiscreteFunctionSpace_);

    const auto& geometry = host_entity.geometry();

    const auto numQuadraturePoints = sg_quadrature.nop();
    for (size_t quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint) {
      const auto& local_point = sg_quadrature.point(quadraturePoint);

      const double quad_weight = sg_quadrature.weight(quadraturePoint) * sg_geometry.integrationElement(local_point);

      // global point in the subgrid
      const auto global_point = sg_geometry.global(local_point);

      // check for subgrid-hostgrid-compatibility
      assert(global_point == geometry.global(quadrature.point(quadraturePoint)));

      std::vector<RangeType> value_b, value_coarse_basis_func, clement_weight;
      value_b.resize(lm_system_matrix.M());
      value_coarse_basis_func.resize(lm_system_matrix.N());
      clement_weight.resize(lm_system_matrix.N());
      for (size_t j = 0; j != lm_system_matrix.M(); ++j) // rows
        ((b_h[j])->localFunction(localgrid_entity)).evaluate(sg_quadrature[quadraturePoint], value_b[j]);
      for (size_t i = 0; i != lm_system_matrix.N(); ++i) // rows
      {
        ((*coarse_basis_)[(*ids_relevant_basis_functions_for_subgrid_)[coarse_index][i]]->localFunction(host_entity))
            .evaluate(quadrature[quadraturePoint], value_coarse_basis_func[i]);
        clement_weight[i] =
            (*inverse_of_L1_norm_coarse_basis_funcs_)[(*ids_relevant_basis_functions_for_subgrid_)[coarse_index][i]];
      }

      for (size_t i = 0; i != lm_system_matrix.N(); ++i)   // rows
        for (size_t j = 0; j != lm_system_matrix.M(); ++j) // colums
          lm_system_matrix[i][j] += quad_weight * clement_weight[i] * value_b[j] * value_coarse_basis_func[i];
    }
  } // lagrange multplier problem system matrix assembled
}

// solve local problem for Local Orthogonal Decomposition Method (LOD)
void LocalProblemSolver::solve_corrector_problem_lod(
    JacobianRangeType& e,
    // the matrix in our linear system of equations
    // in the non-linear case, it is the matrix for each iteration step
    LocProbLinearOperatorTypeType& locprob_system_matrix, MatrixType& lm_system_matrix,
    LocalGridDiscreteFunctionType& local_corrector, const int coarse_index /*= -1*/) const {
  //! if we do not sort out the coarse boundaries (on rigorous_msfem_solver.cc, line 175), the results get better

  // set solution equal to zero:
  local_corrector.clear();

  if (coarse_index < 0)
    DUNE_THROW(Dune::InvalidStateException, "Invalid coarse index: coarse_index < 0");

  const auto& subDiscreteFunctionSpace = local_corrector.space();

  const auto number_of_relevant_coarse_nodes_for_subgrid =
      (*ids_relevant_basis_functions_for_subgrid_)[coarse_index].size();
  const auto& subGrid = subDiscreteFunctionSpace.grid();

  //! define the discrete (elliptic) local MsFEM problem operator
  // ( effect of the discretized differential operator on a certain discrete function )
  LocalProblemOperator local_problem_op(subDiscreteFunctionSpace, diffusion_);

  const auto locprob_inverse_system_matrix = make_inverse_operator(locprob_system_matrix);

  //! Solve the local corrector problem without constraint
  // (without condition that the Lagrange interpolation must be zero)
  // ----------------------------------------------------------------------------------------------------

  // right hand side vectors of the algebraic local MsFEM problem
  LocalGridDiscreteFunctionType corrector_problem_rhs("rhs of local corrector problem in LOD", subDiscreteFunctionSpace);
  corrector_problem_rhs.clear();

  // consider to make separate implementation of 'assemble_local_RHS' for the LOD method
  local_problem_op.assemble_local_RHS(
      e, subgrid_list_.getCoarseNodeVector(coarse_index), /*coarse node vector is a dummy in this case*/
      DSC_CONFIG_GET("msfem.oversampling_strategy", 1),               /*always '3' in this case */
      corrector_problem_rhs);

  local_problem_op.set_zero_boundary_condition_RHS(hostDiscreteFunctionSpace_, corrector_problem_rhs);
  locprob_inverse_system_matrix->apply(corrector_problem_rhs, local_corrector);
  // ----------------------------------------------------------------------------------------------------

  // right hand side vectors for the lagrange multiplier (lm) problems (for e)
  // entries lm_rhs[i] = weight_i ( local_corrector, coarse_basis_func[i] )_L2(\Omega)
  VectorType lm_rhs(number_of_relevant_coarse_nodes_for_subgrid);
  for (size_t i = 0; i != number_of_relevant_coarse_nodes_for_subgrid; ++i) // columns
  {
    lm_rhs[i] = 0.0;
  }

  for (const auto& localgrid_entity : subDiscreteFunctionSpace) {
    const auto& sg_geometry = localgrid_entity.geometry();

    auto host_entity_pointer = subGrid.getLocalEntity<0>(localgrid_entity);
    const auto& host_entity = *host_entity_pointer;

    // exact for polynomials of degree 2:
    const auto sg_quadrature = DSFe::make_quadrature(localgrid_entity, subDiscreteFunctionSpace);
    const auto quadrature = DSFe::make_quadrature(host_entity, hostDiscreteFunctionSpace_);

    RangeType value_local_problem_solution, value_coarse_basis_func_i;

    const auto numQuadraturePoints = sg_quadrature.nop();
    for (size_t quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint) {
      const auto& local_point = sg_quadrature.point(quadraturePoint);

      const double quad_weight = sg_quadrature.weight(quadraturePoint) * sg_geometry.integrationElement(local_point);

      // check for subgrid-hostgrid-compatibility
      assert(sg_geometry.global(local_point) == host_entity.geometry().global(quadrature.point(quadraturePoint)));

      const auto local_sol = local_corrector.localFunction(localgrid_entity);

      local_sol.evaluate(sg_quadrature[quadraturePoint], value_local_problem_solution);

      for (size_t i = 0; i != number_of_relevant_coarse_nodes_for_subgrid; ++i) // columns
      {

        const auto interior_coarse_basis_id_in_subgrid = (*ids_relevant_basis_functions_for_subgrid_)[coarse_index][i];

        const auto local_coarse_basis_i =
            (*coarse_basis_)[interior_coarse_basis_id_in_subgrid]->localFunction(host_entity);
        local_coarse_basis_i.evaluate(quadrature[quadraturePoint], value_coarse_basis_func_i);

        double clement_weight = (*inverse_of_L1_norm_coarse_basis_funcs_)[interior_coarse_basis_id_in_subgrid];

        lm_rhs[i] += quad_weight * clement_weight * value_local_problem_solution * value_coarse_basis_func_i;
      }
    }
  } // lagrange multplier problem system matrix assembled

  MatrixOperatorType lm_matrix_op(lm_system_matrix);

  PreconditionerType lm_preconditioner(lm_system_matrix, 100, 0.9);

  Dune::InverseOperatorResult result_data;
  VectorType v_h(number_of_relevant_coarse_nodes_for_subgrid);
  for (size_t col = 0; col != v_h.N(); ++col) {
    v_h[col] = 0.0;
  }

  typedef Dune::BiCGSTABSolver<VectorType> SolverType;

  double tol = DSC_CONFIG_GET("rigorous_msfem.local_micro_solver_tolerance", 1e-10);
  int num_iterations = DSC_CONFIG_GET("rigorous_msfem.local_micro_solver_iterations", 10000);

  SolverType lm_prob_solver(lm_matrix_op, lm_preconditioner, tol, num_iterations, false);

  lm_prob_solver.apply(v_h, lm_rhs, result_data);

  // set final_rhs = \sum_j clement_weight_j v_h[j] coarse_basis_function[j]
  LocalGridDiscreteFunctionType final_rhs("final rhs local problem", hostDiscreteFunctionSpace_); // for e
  final_rhs.clear();

  for (size_t i = 0; i != number_of_relevant_coarse_nodes_for_subgrid; ++i) // columns
  {

    const auto interior_coarse_basis_id_in_subgrid = (*ids_relevant_basis_functions_for_subgrid_)[coarse_index][i];

    LocalGridDiscreteFunctionType aux_func("auxilliary func", hostDiscreteFunctionSpace_);
    aux_func.clear();

    double coefficient = (*inverse_of_L1_norm_coarse_basis_funcs_)[interior_coarse_basis_id_in_subgrid] * v_h[i];

    aux_func += (*(*coarse_basis_)[interior_coarse_basis_id_in_subgrid]);
    aux_func *= coefficient;
    final_rhs += aux_func;
  }

  // right hand side vectors of the algebraic local MsFEM problem
  LocalGridDiscreteFunctionType final_rhs_vector("final rhs of local MsFEM problem", subDiscreteFunctionSpace); // for e
  final_rhs_vector.clear();

  local_problem_op.assemble_local_RHS_lg_problems(final_rhs, 1.0, final_rhs_vector);
  local_problem_op.set_zero_boundary_condition_RHS(hostDiscreteFunctionSpace_, final_rhs_vector);

  LocalGridDiscreteFunctionType preliminary_solution("preliminary_solution", subDiscreteFunctionSpace); // for e
  preliminary_solution.clear();

  locprob_inverse_system_matrix->apply(final_rhs_vector, preliminary_solution);

  local_corrector -= preliminary_solution;

} // solve_corrector_problem_lod

// solve Dirichlet boundary corrector problem for Local Orthogonal Decomposition Method (LOD)
void LocalProblemSolver::solve_dirichlet_corrector_problem_lod(
    // the matrix in our linear system of equations
    // in the non-linear case, it is the matrix for each iteration step
    LocProbLinearOperatorTypeType& locprob_system_matrix, MatrixType& lm_system_matrix,
    LocalGridDiscreteFunctionType& local_corrector, const int coarse_index /*= -1*/) const {

  // set solution equal to zero:
  local_corrector.clear();

  if (coarse_index < 0)
    DUNE_THROW(Dune::InvalidStateException, "Invalid coarse index: coarse_index < 0");

  const LocalGridDiscreteFunctionSpaceType& subDiscreteFunctionSpace = local_corrector.space();

  auto number_of_relevant_coarse_nodes_for_subgrid = (*ids_relevant_basis_functions_for_subgrid_)[coarse_index].size();
  const LocalGridType& subGrid = subDiscreteFunctionSpace.grid();

  //! define the discrete (elliptic) local MsFEM problem operator
  // ( effect of the discretized differential operator on a certain discrete function )
  LocalProblemOperator local_problem_op(subDiscreteFunctionSpace, diffusion_);

  const auto locprob_inverse_system_matrix = make_inverse_operator(locprob_system_matrix);

  //! Solve the local corrector problem without constraint
  // (without condition that the Lagrange interpolation must be zero)
  // ----------------------------------------------------------------------------------------------------

  // right hand side vectors of the algebraic local MsFEM problem
  LocalGridDiscreteFunctionType corrector_problem_rhs("rhs of local corrector problem in LOD", subDiscreteFunctionSpace);
  corrector_problem_rhs.clear();

  // consider to make separate implementation of 'assemble_local_RHS' for the LOD method
  local_problem_op.assemble_local_RHS_Dirichlet_corrector(
      *dirichlet_extension_,
      subgrid_list_.getCoarseNodeVector(coarse_index), /*coarse node vector is a dummy in this case*/
      DSC_CONFIG_GET("msfem.oversampling_strategy", 1),            /*always '3' in this case */
      corrector_problem_rhs);

  local_problem_op.set_zero_boundary_condition_RHS(hostDiscreteFunctionSpace_, corrector_problem_rhs);
  locprob_inverse_system_matrix->apply(corrector_problem_rhs, local_corrector);

  // ----------------------------------------------------------------------------------------------------

  // right hand side vectors for the lagrange multiplier (lm) problems
  // entries lm_rhs[i] = weight_i ( local_corrector, coarse_basis_func[i] )_L2(\Omega)
  VectorType lm_rhs(number_of_relevant_coarse_nodes_for_subgrid);
  for (size_t i = 0; i != number_of_relevant_coarse_nodes_for_subgrid; ++i) // columns
  {
    lm_rhs[i] = 0.0;
  }

  for (const auto& localgrid_entity : subDiscreteFunctionSpace) {
    const auto& sg_geometry = localgrid_entity.geometry();

    auto host_entity_pointer = subGrid.getLocalEntity<0>(localgrid_entity);
    const auto& host_entity = *host_entity_pointer;

    // exact for polynomials of degree 2:
    const auto sg_quadrature = DSFe::make_quadrature(localgrid_entity, subDiscreteFunctionSpace);
    const auto quadrature = DSFe::make_quadrature(host_entity, hostDiscreteFunctionSpace_);

    RangeType value_local_problem_solution, value_coarse_basis_func_i;

    const auto numQuadraturePoints = sg_quadrature.nop();
    for (size_t quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint) {
      const auto& local_point = sg_quadrature.point(quadraturePoint);

      const double quad_weight = sg_quadrature.weight(quadraturePoint) * sg_geometry.integrationElement(local_point);

      // check for subgrid-hostgrid-compatibility
      assert(sg_geometry.global(local_point) == host_entity.geometry().global(quadrature.point(quadraturePoint)));

      const auto local_sol = local_corrector.localFunction(localgrid_entity);

      local_sol.evaluate(sg_quadrature[quadraturePoint], value_local_problem_solution);

      for (size_t i = 0; i != number_of_relevant_coarse_nodes_for_subgrid; ++i) // columns
      {

        const auto interior_coarse_basis_id_in_subgrid = (*ids_relevant_basis_functions_for_subgrid_)[coarse_index][i];

        auto local_coarse_basis_i = (*coarse_basis_)[interior_coarse_basis_id_in_subgrid]->localFunction(host_entity);
        local_coarse_basis_i.evaluate(quadrature[quadraturePoint], value_coarse_basis_func_i);

        double clement_weight = (*inverse_of_L1_norm_coarse_basis_funcs_)[interior_coarse_basis_id_in_subgrid];

        lm_rhs[i] += quad_weight * clement_weight * value_local_problem_solution * value_coarse_basis_func_i;
      }
    }
  } // lagrange multplier problem system matrix assembled

  MatrixOperatorType lm_matrix_op(lm_system_matrix);

  PreconditionerType lm_preconditioner(lm_system_matrix, 100, 0.9);

  Dune::InverseOperatorResult result_data;
  VectorType v_h(number_of_relevant_coarse_nodes_for_subgrid);
  for (size_t col = 0; col != v_h.N(); ++col) {
    v_h[col] = 0.0;
  }

  typedef Dune::BiCGSTABSolver<VectorType> SolverType;

  double tol = DSC_CONFIG_GET("rigorous_msfem.local_micro_solver_tolerance", 1e-10);
  int num_iterations = DSC_CONFIG_GET("rigorous_msfem.local_micro_solver_iterations", 10000);

  SolverType lm_prob_solver(lm_matrix_op, lm_preconditioner, tol, num_iterations, false);

  lm_prob_solver.apply(v_h, lm_rhs, result_data);

  // set final_rhs = \sum_j clement_weight_j v_h[j] coarse_basis_function[j]
  LocalGridDiscreteFunctionType final_rhs("final rhs local problem", hostDiscreteFunctionSpace_); // for e
  final_rhs.clear();

  for (size_t i = 0; i != number_of_relevant_coarse_nodes_for_subgrid; ++i) // columns
  {
    const auto interior_coarse_basis_id_in_subgrid = (*ids_relevant_basis_functions_for_subgrid_)[coarse_index][i];

    LocalGridDiscreteFunctionType aux_func("auxilliary func", hostDiscreteFunctionSpace_);
    aux_func.clear();

    double coefficient = (*inverse_of_L1_norm_coarse_basis_funcs_)[interior_coarse_basis_id_in_subgrid] * v_h[i];

    aux_func += (*(*coarse_basis_)[interior_coarse_basis_id_in_subgrid]);
    aux_func *= coefficient;
    final_rhs += aux_func;
  }

  // right hand side vectors of the algebraic local MsFEM problem
  LocalGridDiscreteFunctionType final_rhs_vector("final rhs of local MsFEM problem", subDiscreteFunctionSpace); // for e
  final_rhs_vector.clear();

  local_problem_op.assemble_local_RHS_lg_problems(final_rhs, 1.0, final_rhs_vector);
  local_problem_op.set_zero_boundary_condition_RHS(hostDiscreteFunctionSpace_, final_rhs_vector);

  LocalGridDiscreteFunctionType preliminary_solution("preliminary_solution", subDiscreteFunctionSpace); // for e
  preliminary_solution.clear();

  locprob_inverse_system_matrix->apply(final_rhs_vector, preliminary_solution);

  local_corrector -= preliminary_solution;

} // solve_dirichlet_corrector_problem_lod

// solve Neumann boundary corrector problem for Local Orthogonal Decomposition Method (LOD)
void LocalProblemSolver::solve_neumann_corrector_problem_lod(
    // the matrix in our linear system of equations
    // in the non-linear case, it is the matrix for each iteration step
    LocProbLinearOperatorTypeType& locprob_system_matrix, MatrixType& lm_system_matrix,
    LocalGridDiscreteFunctionType& local_corrector, const int coarse_index /*= -1*/) const {

  // set solution equal to zero:
  local_corrector.clear();

  if (coarse_index < 0)
    DUNE_THROW(Dune::InvalidStateException, "Invalid coarse index: coarse_index < 0");

  const auto& subDiscreteFunctionSpace = local_corrector.space();

  const auto number_of_relevant_coarse_nodes_for_subgrid = (*ids_relevant_basis_functions_for_subgrid_)[coarse_index].size();
  const auto& subGrid = subDiscreteFunctionSpace.grid();

  //! define the discrete (elliptic) local MsFEM problem operator
  // ( effect of the discretized differential operator on a certain discrete function )
  LocalProblemOperator local_problem_op(subDiscreteFunctionSpace, diffusion_);

  const auto locprob_inverse_system_matrix = make_inverse_operator(locprob_system_matrix);

  //! Solve the local corrector problem without constraint
  // (without condition that the Lagrange interpolation must be zero)
  // ----------------------------------------------------------------------------------------------------

  // right hand side vectors of the algebraic local MsFEM problem
  LocalGridDiscreteFunctionType corrector_problem_rhs("rhs of local corrector problem in LOD", subDiscreteFunctionSpace);
  corrector_problem_rhs.clear();

  // consider to make separate implementation of 'assemble_local_RHS' for the LOD method
  local_problem_op.assemble_local_RHS_Neumann_corrector(
      *neumann_bc_, hostDiscreteFunctionSpace_,
      subgrid_list_.getCoarseNodeVector(coarse_index), /*coarse node vector is a dummy in this case*/
      DSC_CONFIG_GET("msfem.oversampling_strategy", 1),            /*always '3' in this case */
      corrector_problem_rhs);

  local_problem_op.set_zero_boundary_condition_RHS(hostDiscreteFunctionSpace_, corrector_problem_rhs);
  locprob_inverse_system_matrix->apply(corrector_problem_rhs, local_corrector);
  // ----------------------------------------------------------------------------------------------------

  // right hand side vectors for the lagrange multiplier (lm) problems
  // entries lm_rhs[i] = weight_i ( local_corrector, coarse_basis_func[i] )_L2(\Omega)
  VectorType lm_rhs(number_of_relevant_coarse_nodes_for_subgrid);
  for (size_t i = 0; i != number_of_relevant_coarse_nodes_for_subgrid; ++i) // columns
  {
    lm_rhs[i] = 0.0;
  }

  for (const auto& localgrid_entity : subDiscreteFunctionSpace) {
    const auto& sg_geometry = localgrid_entity.geometry();

    const auto host_entity_pointer = subGrid.getLocalEntity<0>(localgrid_entity);
    const auto& host_entity = *host_entity_pointer;

    // exact for polynomials of degree 2:
    const auto sg_quadrature = DSFe::make_quadrature(localgrid_entity, subDiscreteFunctionSpace);
    const auto quadrature = DSFe::make_quadrature(host_entity, hostDiscreteFunctionSpace_);

    RangeType value_local_problem_solution, value_coarse_basis_func_i;

    const auto numQuadraturePoints = sg_quadrature.nop();
    for (size_t quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint) {
      const auto& local_point = sg_quadrature.point(quadraturePoint);

      const double quad_weight = sg_quadrature.weight(quadraturePoint) * sg_geometry.integrationElement(local_point);

      // check for subgrid-hostgrid-compatibility
      assert(sg_geometry.global(local_point) == host_entity.geometry().global(quadrature.point(quadraturePoint)));

      const auto local_sol = local_corrector.localFunction(localgrid_entity);

      local_sol.evaluate(sg_quadrature[quadraturePoint], value_local_problem_solution);

      for (size_t i = 0; i != number_of_relevant_coarse_nodes_for_subgrid; ++i) // columns
      {

        const auto interior_coarse_basis_id_in_subgrid = (*ids_relevant_basis_functions_for_subgrid_)[coarse_index][i];

        const auto local_coarse_basis_i =
            (*coarse_basis_)[interior_coarse_basis_id_in_subgrid]->localFunction(host_entity);
        local_coarse_basis_i.evaluate(quadrature[quadraturePoint], value_coarse_basis_func_i);

        const double clement_weight = (*inverse_of_L1_norm_coarse_basis_funcs_)[interior_coarse_basis_id_in_subgrid];

        lm_rhs[i] += quad_weight * clement_weight * value_local_problem_solution * value_coarse_basis_func_i;
      }
    }
  } // lagrange multplier problem system matrix assembled

  MatrixOperatorType lm_matrix_op(lm_system_matrix);
  PreconditionerType lm_preconditioner(lm_system_matrix, 100, 0.9);
  Dune::InverseOperatorResult result_data;

  VectorType v_h(number_of_relevant_coarse_nodes_for_subgrid);
  for (size_t col = 0; col != v_h.N(); ++col) {
    v_h[col] = 0.0;
  }

  typedef Dune::BiCGSTABSolver<VectorType> SolverType;

  const double tol = DSC_CONFIG_GET("rigorous_msfem.local_micro_solver_tolerance", 1e-10);
  const int num_iterations = DSC_CONFIG_GET("rigorous_msfem.local_micro_solver_iterations", 10000);

  SolverType lm_prob_solver(lm_matrix_op, lm_preconditioner, tol, num_iterations, false);

  lm_prob_solver.apply(v_h, lm_rhs, result_data);

  // set final_rhs = \sum_j clement_weight_j v_h[j] coarse_basis_function[j]
  LocalGridDiscreteFunctionType final_rhs("final rhs local problem", hostDiscreteFunctionSpace_); // for e
  final_rhs.clear();

  for (size_t i = 0; i != number_of_relevant_coarse_nodes_for_subgrid; ++i) // columns
  {

    const auto interior_coarse_basis_id_in_subgrid = (*ids_relevant_basis_functions_for_subgrid_)[coarse_index][i];

    LocalGridDiscreteFunctionType aux_func("auxilliary func", hostDiscreteFunctionSpace_);
    aux_func.clear();

    const double coefficient = (*inverse_of_L1_norm_coarse_basis_funcs_)[interior_coarse_basis_id_in_subgrid] * v_h[i];

    aux_func += (*(*coarse_basis_)[interior_coarse_basis_id_in_subgrid]);
    aux_func *= coefficient;
    final_rhs += aux_func;
  }

  // right hand side vectors of the algebraic local MsFEM problem
  LocalGridDiscreteFunctionType final_rhs_vector("final rhs of local MsFEM problem", subDiscreteFunctionSpace); // for e
  final_rhs_vector.clear();

  local_problem_op.assemble_local_RHS_lg_problems(final_rhs, 1.0, final_rhs_vector);
  local_problem_op.set_zero_boundary_condition_RHS(hostDiscreteFunctionSpace_, final_rhs_vector);

  LocalGridDiscreteFunctionType preliminary_solution("preliminary_solution", subDiscreteFunctionSpace); // for e
  preliminary_solution.clear();

  locprob_inverse_system_matrix->apply(final_rhs_vector, preliminary_solution);

  local_corrector -= preliminary_solution;
} // solve_dirichlet_neumann_problem_lod
#endif // ENABLE_LOD_ONLY_CODE

} // namespace MsFEM {
} // namespace Multiscale {
} // namespace Dune {
