#include <config.h>
#include <assert.h>
#include <boost/concept/usage.hpp>
#include <boost/format.hpp>
#include <dune/common/exceptions.hh>
#include <dune/istl/scalarproducts.hh>
#include <dune/istl/solvers.hh>
#include <dune/multiscale/msfem/localproblems/localoperator.hh>
#include <dune/multiscale/problems/selector.hh>
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
#include <dune/stuff/functions/norm.hh>
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

/** \brief define output parameters for local problems
 *  appends "local_problems" for path
 **/
struct LocalProblemDataOutputParameters : public OutputParameters {
public:
  explicit LocalProblemDataOutputParameters();
};

LocalProblemDataOutputParameters::LocalProblemDataOutputParameters()
  : OutputParameters(DSC_CONFIG_GET("global.datadir", "data") + "/local_problems/") {}

LocalProblemSolver::LocalProblemSolver(const CommonTraits::DiscreteFunctionSpaceType& coarse_space,
                                       LocalGridList& subgrid_list,
                                       const CommonTraits::DiffusionType& diffusion_operator)
  : diffusion_(diffusion_operator)
  , subgrid_list_(subgrid_list)
  , coarse_space_(coarse_space) {}


void LocalProblemSolver::solve_all_on_single_cell(const MsFEMTraits::CoarseEntityType& coarseCell,
                                                  MsFEMTraits::LocalSolutionVectorType& allLocalSolutions) const {
  assert(allLocalSolutions.size() > 0);

  const bool hasBoundary = coarseCell.hasBoundaryIntersections();
  const auto numBoundaryCorrectors = DSG::is_simplex_grid(coarse_space_) ? 1u : 2u;
  const auto numInnerCorrectors = allLocalSolutions.size() - numBoundaryCorrectors;

  // clear return argument
  for (auto& localSol : allLocalSolutions)
    localSol->vector() *= 0;

  const auto& subDiscreteFunctionSpace = allLocalSolutions[0]->space();

  //! define the discrete (elliptic) local MsFEM problem operator
  // ( effect of the discretized differential operator on a certain discrete function )
  LocalProblemOperator localProblemOperator(coarse_space_, subDiscreteFunctionSpace, diffusion_);

  // right hand side vector of the algebraic local MsFEM problem
  MsFEMTraits::LocalSolutionVectorType allLocalRHS(allLocalSolutions.size());
  for (auto& it : allLocalRHS)
    it = DSC::make_unique<MsFEMTraits::LocalGridDiscreteFunctionType>(subDiscreteFunctionSpace,
                                                                      "rhs of local MsFEM problem");

  localProblemOperator.assemble_all_local_rhs(coarseCell, allLocalRHS);

  for (auto i : DSC::valueRange(allLocalSolutions.size())) {
    auto& current_rhs = *allLocalRHS[i];
    auto& current_solution = *allLocalSolutions[i];

    // is the right hand side of the local MsFEM problem equal to zero or almost identical to zero?
    // if yes, the solution of the local MsFEM problem is also identical to zero. The solver is getting a problem with
    // this situation, which is why we do not solve local msfem problems for zero-right-hand-side, since we already know
    // the result.
    const auto norm = GDT::Products::L2< typename MsFEMTraits::LocalGridViewType >
                      (*current_rhs.space().grid_view()).induced_norm(current_rhs);
    if (norm < 1e-12) {
      current_solution.vector() *= 0;
      DSC_LOG_DEBUG << "Local MsFEM problem with solution zero." << std::endl;
      continue;
    }
    // don't solve local problems for boundary correctors if coarse cell has no boundary intersections
    if (i >= numInnerCorrectors && !hasBoundary) {
      current_solution.vector() *= 0;
      DSC_LOG_DEBUG << "Zero-Boundary corrector." << std::endl;
      continue;
    }
    localProblemOperator.apply_inverse(current_rhs, current_solution);
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
const auto& grid = coarse_space_.grid_view()->grid();
  const auto coarseGridSize = grid.size(0) - grid.overlapSize(0);
  if (Dune::Fem::MPIManager::size() > 0)
    DSC_LOG_DEBUG << "Rank " << Dune::Fem::MPIManager::rank() << " will solve local problems for " << coarseGridSize
                 << " coarse entities!" << std::endl;
  else {
    DSC_LOG_DEBUG << "Will solve local problems for " << coarseGridSize << " coarse entities!" << std::endl;
  }
  DSC_PROFILER.startTiming("msfem.localProblemSolver.solveAndSaveAll");

  // we want to determine minimum, average and maxiumum time for solving a local msfem problem in the current method
  DSC::MinMaxAvg<double> solveTime;

  const auto& coarseGridLeafIndexSet = coarse_space_.grid_view()->grid().leafIndexSet();
  for (const auto& coarseEntity : DSC::viewRange(*coarse_space_.grid_view())) {
    const int coarse_index = coarseGridLeafIndexSet.index(coarseEntity);
    DSC_LOG_DEBUG << "-------------------------" << std::endl << "Coarse index " << coarse_index << std::endl;

    // take time
    DSC_PROFILER.startTiming("msfem.localProblemSolver.solve");
    LocalSolutionManager localSolutionManager(coarse_space_, coarseEntity, subgrid_list_);
    // solve the problems
    solve_all_on_single_cell(coarseEntity, localSolutionManager.getLocalSolutions());
    solveTime(DSC_PROFILER.stopTiming("msfem.localProblemSolver.solve") / 1000.f);

    // save the local solutions to disk/mem
    localSolutionManager.save();

    DSC_PROFILER.resetTiming("msfem.localProblemSolver.solve");
  } // for

  //! @todo The following debug-output is wrong (number of local problems may be different)
  const auto totalTime = DSC_PROFILER.stopTiming("msfem.localProblemSolver.solveAndSaveAll") / 1000.f;
  DSC_LOG_DEBUG << "Local problems solved for " << coarseGridSize << " coarse grid entities.\n";
  DSC_LOG_DEBUG << "Minimum time for solving a local problem = " << solveTime.min() << "s.\n";
  DSC_LOG_DEBUG << "Maximum time for solving a local problem = " << solveTime.max() << "s.\n";
  DSC_LOG_DEBUG << "Average time for solving a local problem = " << solveTime.average() << "s.\n";
  DSC_LOG_DEBUG << "Total time for computing and saving the localproblems = "
                << totalTime << "s on rank" << Dune::Fem::MPIManager::rank() << std::endl;
} // assemble_all

} // namespace MsFEM {
} // namespace Multiscale {
} // namespace Dune {
