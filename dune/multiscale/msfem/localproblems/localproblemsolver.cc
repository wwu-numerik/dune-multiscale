#include <config.h>
#include <assert.h>
#include <boost/concept/usage.hpp>
#include <boost/format.hpp>
#include <dune/common/exceptions.hh>
#include <dune/istl/scalarproducts.hh>
#include <dune/istl/solvers.hh>
#include <dune/multiscale/msfem/localproblems/localoperator.hh>
#include <dune/multiscale/problems/selector.hh>
#include <dune/multiscale/common/df_io.hh>
#include <dune/multiscale/msfem/localproblems/localsolutionmanager.hh>
#include <dune/multiscale/tools/misc.hh>
#include <dune/xt/common/logging.hh>
#include <dune/xt/common/math.hh>
#include <dune/xt/common/memory.hh>
#include <dune/xt/common/configuration.hh>
#include <dune/xt/common/timings.hh>
#include <dune/xt/common/ranges.hh>
#include <dune/xt/common/parallel/partitioner.hh>
#include <dune/xt/grid/parallel/partitioning/ranged.hh>
#include <dune/gdt/operators/l2.hh>
#include <dune/xt/grid/walker.hh>
#include <dune/xt/grid/walker/functors.hh>
#include <iterator>
#include <memory>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include <dune/multiscale/msfem/localproblems/localgridlist.hh>
#include <dune/multiscale/tools/misc/outputparameter.hh>
#include "localproblemsolver.hh"

namespace Dune {
namespace Multiscale {

/** \brief define output parameters for local problems
 *  appends "local_problems" for path
 **/
struct LocalProblemDataOutputParameters : public OutputParameters
{
public:
  explicit LocalProblemDataOutputParameters(const Problem::ProblemContainer& problem);
};

LocalProblemDataOutputParameters::LocalProblemDataOutputParameters(const DMP::ProblemContainer& problem)
  : OutputParameters(problem.config().get("global.datadir", "data") + std::string("/local_problems/"))
{
}

LocalProblemSolver::LocalProblemSolver(const Problem::ProblemContainer& problem,
                                       CommonTraits::SpaceType coarse_space,
                                       LocalGridList& localgrid_list)
  : localgrid_list_(localgrid_list)
  , coarse_space_(coarse_space)
  , problem_(problem)
{
}

void LocalProblemSolver::solve_all_on_single_cell(
    const MsFEMTraits::CoarseEntityType& coarseCell,
    MsFEMTraits::LocalSolutionVectorType& all_localproblem_solutions) const
{
  assert(all_localproblem_solutions.size() > 0);

  const bool hasBoundary = coarseCell.hasBoundaryIntersections();
  const auto numBoundaryCorrectors = Dune::XT::Grid::is_simplex_grid(*coarse_space_) ? 1u : 2u;
  const auto numInnerCorrectors = all_localproblem_solutions.size() - numBoundaryCorrectors;

  // clear return argument
  for (auto& localSol : all_localproblem_solutions)
    localSol->vector() *= 0;

  const auto& local_space = all_localproblem_solutions[0]->space();

  //! define the discrete (elliptic) local MsFEM problem operator
  // ( effect of the discretized differential operator on a certain discrete function )
  LocalProblemOperator localProblemOperator(problem_, *coarse_space_, local_space);

  // right hand side vector of the algebraic local MsFEM problem
  MsFEMTraits::LocalSolutionVectorType allLocalRHS(all_localproblem_solutions.size());
  for (auto& it : allLocalRHS)
    it = Dune::XT::Common::make_unique<MsFEMTraits::LocalGridDiscreteFunctionType>(local_space,
                                                                                   "rhs of local MsFEM problem");

  localProblemOperator.assemble_all_local_rhs(coarseCell, allLocalRHS);

  for (auto i : Dune::XT::Common::value_range(all_localproblem_solutions.size())) {
    auto& current_rhs = *allLocalRHS[i];
    auto& current_solution = *all_localproblem_solutions[i];

    // is the right hand side of the local MsFEM problem equal to zero or almost identical to zero?
    // if yes, the solution of the local MsFEM problem is also identical to zero. The solver is getting a problem with
    // this situation, which is why we do not solve local msfem problems for zero-right-hand-side, since we already know
    // the result.
    //! TODO calculating the norm has really bad performance impact, is the instability actually still there?
    //    const auto norm = GDT::Products::L2<typename MsFEMTraits::LocalGridViewType>(current_rhs.space().grid_view())
    //                          .induced_norm(current_rhs);
    //    if (norm < 1e-12) {
    //      current_solution.vector() *= 0;
    //      MS_LOG_DEBUG << boost::format("Local MsFEM problem with solution zero. (corrector %d)") % i << std::endl;
    //      continue;
    //    }

    // don't solve local problems for boundary correctors if coarse cell has no boundary intersections
    if (i >= numInnerCorrectors && !hasBoundary) {
      current_solution.vector() *= 0;
      MS_LOG_DEBUG << "Zero-Boundary corrector." << std::endl;
      continue;
    }
    localProblemOperator.apply_inverse(current_rhs, current_solution);
  }
}

void LocalProblemSolver::solve_for_all_cells()
{
  const auto& grid = coarse_space_->grid_layer().grid();
  const auto coarseGridSize = grid.size(0) - grid.overlapSize(0);

  MS_LOG_INFO << boost::format("Rank %d will solve local problems for %d coarse entities\n") % grid.comm().rank()
                     % coarseGridSize;
  DXTC_TIMINGS.start("msfem.local.solve_for_all_cells");

  // we want to determine minimum, average and maxiumum time for solving a local msfem problem in the current method
  Dune::XT::Common::MinMaxAvg<double> solveTime;

  const auto interior = coarse_space_->grid_layer().grid().leafGridView();
  typedef std::remove_const<decltype(interior)>::type InteriorType;
  GDT::SystemAssembler<CommonTraits::SpaceType, InteriorType> walker(*coarse_space_, interior);

  const auto threads = Dune::XT::Common::threadManager().max_threads();
  XT::Grid::RangedPartitioning<InteriorType, 0, CommonTraits::InteriorBorderPartition> partitioning(interior, threads);

  const std::function<void(const CommonTraits::EntityType&)> func = [&](const CommonTraits::EntityType& coarseEntity) {
    const int coarse_index = walker.ansatz_space().grid_layer().indexSet().index(coarseEntity);
    MS_LOG_DEBUG << "-------------------------" << std::endl << "Coarse index " << coarse_index << std::endl;

    // take time
    //    DXTC_TIMINGS.start("msfem.local.solve_all_on_single_cell");
    LocalproblemSolutionManager localSolutionManager(*coarse_space_, coarseEntity, localgrid_list_);
    // solve the problems
    solve_all_on_single_cell(coarseEntity, localSolutionManager.getLocalSolutions());
    //    solveTime(DXTC_TIMINGS.stop("msfem.local.solve_all_on_single_cell") / 1000.f);

    // save the local solutions to disk/mem
    localSolutionManager.save();

    //    DXTC_TIMINGS.resetTiming("msfem.local.solve_all_on_single_cell");
  };

  walker.add(func);
  walker.assemble(partitioning);

  //! @todo The following debug-output is wrong (number of local problems may be different)
  const auto totalTime = DXTC_TIMINGS.stop("msfem.local.solve_for_all_cells") / 1000.f;
  MS_LOG_INFO << "Local problems solved for " << coarseGridSize << " coarse grid entities.\n"
              //               << "Minimum time for solving a local problem = " << solveTime.min() << "s.\n"
              //               << "Maximum time for solving a local problem = " << solveTime.max() << "s.\n"
              //               << "Average time for solving a local problem = " << solveTime.average() << "s.\n"
              << "Total time for computing and saving the localproblems = " << totalTime << "s on rank"
              << coarse_space_->grid_layer().grid().comm().rank() << std::endl;
} // assemble_all

} // namespace Multiscale {
} // namespace Dune {
