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
#include <dune/multiscale/tools/misc.hh>
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

std::unique_ptr<LocalProblemSolver::InverseOperatorType>
LocalProblemSolver::make_inverse_operator(LocalProblemSolver::LinearOperatorTypeType& problem_matrix) {
  const auto solver =
      Dune::Multiscale::Problem::getModelData()->symmetricDiffusion() ? std::string("cg") : std::string("bcgs");
  return DSC::make_unique<InverseOperatorType>(
      problem_matrix, 1e-8, 1e-8, 20000, DSC_CONFIG_GET("msfem.localproblemsolver_verbose", false), solver,
      DSC_CONFIG_GET("preconditioner_type", std::string("sor")), 1);
}

LocalProblemSolver::LocalProblemSolver(const CommonTraits::DiscreteFunctionSpaceType &coarse_space, LocalGridList& subgrid_list,
                                                 const CommonTraits::DiffusionType &diffusion_operator)
  : diffusion_(diffusion_operator)
  , subgrid_list_(subgrid_list)
  , coarse_space_(coarse_space)
{}


/** Solve all local MsFEM problems for one coarse entity at once.
*
*
*/
void LocalProblemSolver::solve_all_on_single_cell(const MsFEMTraits::CoarseEntityType& coarseCell,
                                                    MsFEMTraits::LocalSolutionVectorType &allLocalSolutions) const {
  assert(allLocalSolutions.size() > 0);

  const bool hasBoundary = coarseCell.hasBoundaryIntersections();
  const auto numBoundaryCorrectors = DSG::is_simplex_grid(coarse_space_) ? 1u : 2u;
  const auto numInnerCorrectors = allLocalSolutions.size() - numBoundaryCorrectors;

  // clear return argument
  for (auto& localSol : allLocalSolutions)
    localSol->clear();

  const auto& subDiscreteFunctionSpace = allLocalSolutions[0]->space();

  //! the matrix in our linear system of equations
  // in the non-linear case, it is the matrix for each iteration step
  LinearOperatorTypeType locProbSysMatrix("Local Problem System Matrix", subDiscreteFunctionSpace,
                                                 subDiscreteFunctionSpace);

  //! define the discrete (elliptic) local MsFEM problem operator
  // ( effect of the discretized differential operator on a certain discrete function )
  LocalProblemOperator localProblemOperator(coarse_space_, subDiscreteFunctionSpace, diffusion_);

  // right hand side vector of the algebraic local MsFEM problem
  MsFEMTraits::LocalSolutionVectorType allLocalRHS(allLocalSolutions.size());
  for (auto& it : allLocalRHS)
    it = DSC::make_unique<MsFEMTraits::LocalGridDiscreteFunctionType>("rhs of local MsFEM problem", subDiscreteFunctionSpace);

  localProblemOperator.assemble_matrix(locProbSysMatrix);
  localProblemOperator.assemble_all_local_rhs(coarseCell,allLocalRHS);

  // set dirichlet dofs to zero
  Stuff::GridboundaryAllDirichlet<MsFEMTraits::LocalGridType::LeafGridView::Intersection> boundaryInfo;
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
}


void LocalProblemSolver::solve_for_all_cells() {
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
  const auto coarseGridSize = coarse_space_.grid().size(0);
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

  const auto& coarseGridLeafIndexSet = coarse_space_.gridPart().grid().leafIndexSet();
  for (const auto& coarseEntity : coarse_space_) {
    const int coarse_index = coarseGridLeafIndexSet.index(coarseEntity);

    DSC_LOG_INFO << "-------------------------" << std::endl << "Coarse index " << coarse_index << std::endl;

    DSC_PROFILER.startTiming("none.saveLocalProblemsOnCell");

    // take time
    DSC_PROFILER.startTiming("none.local_problem_solution");
    LocalSolutionManager localSolutionManager(coarse_space_, coarseEntity, subgrid_list_);

    // solve the problems
    solve_all_on_single_cell(coarseEntity, localSolutionManager.getLocalSolutions());
    // min/max time
    cell_time(DSC_PROFILER.stopTiming("none.local_problem_solution") / 1000.f);
    DSC_PROFILER.resetTiming("none.local_problem_solution");

    // save the local solutions to disk
    DSC_PROFILER.startTiming("none.saveLocalProblemSolution");
    localSolutionManager.save();
    saveTime(DSC_PROFILER.stopTiming("none.saveLocalProblemSolution") / 1000.f);
    DSC_PROFILER.resetTiming("none.saveLocalProblemSolution");


    DSC_LOG_INFO << "Total time for solving and saving all local problems for the current subgrid: "
                 << DSC_PROFILER.stopTiming("none.saveLocalProblemsOnCell") / 1000.f << "s"
                 << std::endl << std::endl;
    DSC_PROFILER.resetTiming("none.saveLocalProblemsOnCell");
  } // for

  //! @todo The following debug-output is wrong (number of local problems may be different)
  const auto total_time = DSC_PROFILER.stopTiming("msfem.localproblemsolver.assemble_all") / 1000.f;
  DSC_LOG_INFO << std::endl;
  DSC_LOG_INFO << "Local problems solved for " << coarseGridSize << " coarse grid entities.\n";
  DSC_LOG_INFO << "Minimum time for solving a local problem = " << cell_time.min() << "s.\n";
  DSC_LOG_INFO << "Maximum time for solving a localproblem = " << cell_time.max() << "s.\n";
  DSC_LOG_INFO << "Average time for solving a localproblem = " << cell_time.average() << "s.\n";
  DSC_LOG_INFO << "Minimum time for saving a local problem = " << saveTime.min() << "s.\n";
  DSC_LOG_INFO << "Maximum time for saving a localproblem = " << saveTime.max() << "s.\n";
  DSC_LOG_INFO << "Average time for saving a localproblem = " << saveTime.average() << "s.\n";
  DSC_LOG_INFO << "Total time for computing and saving the localproblems = " << total_time << "s," << std::endl
               << std::endl;
} // assemble_all

} // namespace MsFEM {
} // namespace Multiscale {
} // namespace Dune {
