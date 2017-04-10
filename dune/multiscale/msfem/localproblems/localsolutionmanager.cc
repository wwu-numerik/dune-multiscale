#include <config.h>
#include "localsolutionmanager.hh"

#include <boost/format.hpp>
#include <dune/xt/common/memory.hh>
#include <dune/multiscale/common/df_io.hh>
#include <dune/multiscale/tools/misc.hh>
#include <dune/multiscale/msfem/localproblems/localgridlist.hh>
#include <dune/multiscale/common/df_io.hh>

namespace Dune {
namespace Multiscale {

LocalproblemSolutionManager::LocalproblemSolutionManager(const CommonTraits::SpaceType& coarse_space,
                                                         const MsFEMTraits::CoarseEntityType& coarseEntity,
                                                         const LocalGridList& subgridList)
  : subgridList_(subgridList)
  , subgrid_(subgridList_.getSubGrid(coarseEntity))
  , grid_view_(subgrid_.leafGridView())
  , numBoundaryCorrectors_(DSG::is_simplex_grid(coarse_space) ? 1 : 2)
  , numLocalProblems_(DSG::is_simplex_grid(coarse_space) ? CommonTraits::world_dim + 1
                                                         : coarse_space.mapper().maxNumDofs() + 2)
  , localSolutions_(numLocalProblems_)
  , localSolutionLocation_((boost::format("local_problems/_localProblemSolutions_%d")
                            % coarse_space.grid_view().grid().leafIndexSet().index(coarseEntity))
                               .str())
  , memory_backend_(DiscreteFunctionIO::memory(localSolutionLocation_, grid_view_))
{
  for (auto& it : localSolutions_)
    it = make_df_ptr<MsFEMTraits::LocalGridDiscreteFunctionType>("Local problem Solution", memory_backend_.space());
}

MsFEMTraits::LocalSolutionVectorType& LocalproblemSolutionManager::getLocalSolutions()
{
  return localSolutions_;
}

const MsFEMTraits::LocalSpaceType& LocalproblemSolutionManager::space() const
{
  return memory_backend_.space();
}

void LocalproblemSolutionManager::load()
{
  assert(localSolutions_.size() >= numLocalProblems_);
  for (unsigned int i = 0; i < numLocalProblems_; ++i) {
    auto& solution = localSolutions_[i];
    assert(solution);
    solution->vector() *= 0;
    memory_backend_.read(i, solution);
  }
} // load

void LocalproblemSolutionManager::save() const
{
  for (auto& it : localSolutions_)
    memory_backend_.append(it);
} // save

std::size_t LocalproblemSolutionManager::numBoundaryCorrectors() const
{
  return numBoundaryCorrectors_;
}

} // namespace Multiscale {
} // namespace Dune {
