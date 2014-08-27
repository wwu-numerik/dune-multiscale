#include <config.h>
#include <boost/format.hpp>
#include <dune/stuff/common/memory.hh>
#include <dune/multiscale/tools/discretefunctionwriter.hh>
#include <dune/multiscale/tools/misc.hh>

#include "localsolutionmanager.hh"

namespace Dune {
namespace Multiscale {
namespace MsFEM {

typedef DiscreteFunctionIO<MsFEMTraits::LocalGridDiscreteFunctionType> IOType;

LocalSolutionManager::LocalSolutionManager(const CommonTraits::DiscreteFunctionSpaceType& coarse_space,
                                           const CoarseEntityType& coarseEntity, const LocalGridList &subgridList)
  : subgridList_(subgridList)
  , subgrid_(subgridList_.getSubGrid(coarseEntity))
  , grid_view_ptr_(std::make_shared<MsFEMTraits::LocalGridViewType>(subgrid_.leafGridView()))
  , numBoundaryCorrectors_(DSG::is_simplex_grid(coarse_space) ? 1 : 2)
  , numLocalProblems_(DSG::is_simplex_grid(coarse_space) ? GridSelector::dimgrid + 1
                                                         : coarse_space.mapper().maxNumDofs() + 2)
  , localSolutions_(numLocalProblems_)
  , localSolutionLocation_((boost::format("local_problems/_localProblemSolutions_%d") %
                            coarse_space.grid_view()->grid().leafIndexSet().index(coarseEntity)).str()) {

  auto& reader = IOType::memory(localSolutionLocation_, grid_view_ptr_);
  for (auto& it : localSolutions_)
    it = make_df_ptr<LocalGridDiscreteFunctionType>("Local problem Solution", reader.space());
}

MsFEMTraits::LocalSolutionVectorType& LocalSolutionManager::getLocalSolutions() { return localSolutions_; }

const LocalSolutionManager::LocalGridDiscreteFunctionSpaceType& LocalSolutionManager::space() const {
  return IOType::memory(localSolutionLocation_, grid_view_ptr_).space();
}

void LocalSolutionManager::load() {
  auto& reader = IOType::memory(localSolutionLocation_, grid_view_ptr_);
  for (unsigned int i = 0; i < numLocalProblems_; ++i) {
    localSolutions_[i]->vector() *= 0;
    reader.read(i, localSolutions_[i]);
  }
} // load

void LocalSolutionManager::save() const {
  auto& writer = IOType::memory(localSolutionLocation_, grid_view_ptr_);
  for (auto& it : localSolutions_)
    writer.append(it);
} // save

std::size_t LocalSolutionManager::numBoundaryCorrectors() const { return numBoundaryCorrectors_; }

} // namespace MsFEM {
} // namespace Multiscale {
} // namespace Dune {
