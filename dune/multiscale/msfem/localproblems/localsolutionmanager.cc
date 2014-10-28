#include <config.h>
#include <boost/format.hpp>
#include <dune/stuff/common/memory.hh>
#include <dune/multiscale/tools/discretefunctionwriter.hh>
#include <dune/multiscale/tools/misc.hh>

#include "localsolutionmanager.hh"

namespace Dune {
namespace Multiscale {


LocalSolutionManager::LocalSolutionManager(const CommonTraits::SpaceType& coarse_space,
                                           const CoarseEntityType& coarseEntity, const LocalGridList& subgridList)
  : subgridList_(subgridList)
  , subgrid_(subgridList_.getSubGrid(coarseEntity))
  , grid_view_(subgrid_.leafGridView())
  , numBoundaryCorrectors_(DSG::is_simplex_grid(coarse_space) ? 1 : 2)
  , numLocalProblems_(DSG::is_simplex_grid(coarse_space) ? GridSelector::dimgrid + 1
                                                         : coarse_space.mapper().maxNumDofs() + 2)
  , localSolutions_(numLocalProblems_)
  , localSolutionLocation_((boost::format("local_problems/_localProblemSolutions_%d") %
                            coarse_space.grid_view().grid().leafIndexSet().index(coarseEntity)).str())
  , memory_backend_(IOType::memory(localSolutionLocation_, grid_view_))
{
  for (auto& it : localSolutions_)
    it = make_df_ptr<LocalGridDiscreteFunctionType>("Local problem Solution", memory_backend_.space());
}

MsFEMTraits::LocalSolutionVectorType& LocalSolutionManager::getLocalSolutions() { return localSolutions_; }

const LocalSolutionManager::LocalSpaceType& LocalSolutionManager::space() const {
  return memory_backend_.space();
}

void LocalSolutionManager::load() {
  for (unsigned int i = 0; i < numLocalProblems_; ++i) {
    localSolutions_[i]->vector() *= 0;
    memory_backend_.read(i, localSolutions_[i]);
  }
} // load

void LocalSolutionManager::save() const {
  for (auto& it : localSolutions_)
    memory_backend_.append(it);
} // save

std::size_t LocalSolutionManager::numBoundaryCorrectors() const { return numBoundaryCorrectors_; }

} // namespace Multiscale {
} // namespace Dune {
