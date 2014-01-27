#include <config.h>
#include <boost/format.hpp>
#include <dune/stuff/common/memory.hh>
#include <dune/multiscale/tools/discretefunctionwriter.hh>

#include "localsolutionmanager.hh"

namespace Dune {
namespace Multiscale {
namespace MsFEM {

typedef DiscreteFunctionIO<MsFEMTraits::LocalGridDiscreteFunctionType> IOType;

LocalSolutionManager::LocalSolutionManager(const CoarseEntityType& coarseEntity, LocalGridList& subgridList,
                                           const MacroMicroGridSpecifier& gridSpecifier)
  : subgridList_(subgridList)
  , subgrid_(subgridList_.getSubGrid(coarseEntity))
  , gridSpecifier_(gridSpecifier)
  , numBoundaryCorrectors_((gridSpecifier_.simplexCoarseGrid()) ? 1 : 2)
  , numLocalProblems_((gridSpecifier_.simplexCoarseGrid()) ? GridSelector::dimgrid + 1
                                                           : gridSpecifier_.coarseSpace().mapper().maxNumDofs() + 2)
  , localSolutions_(numLocalProblems_)
  , localSolutionLocation_((boost::format("local_problems/_localProblemSolutions_%d")
                            % gridSpecifier_.coarseSpace().gridPart().grid().leafIndexSet().index(coarseEntity)).str()) {
  auto& reader = IOType::memory(localSolutionLocation_, subgrid_);
  for (auto& it : localSolutions_)
    it = make_df_ptr<LocalGridDiscreteFunctionType>("Local problem Solution", reader.space());
}

MsFEMTraits::LocalSolutionVectorType& LocalSolutionManager::getLocalSolutions() { return localSolutions_; }

const LocalSolutionManager::LocalGridDiscreteFunctionSpaceType& LocalSolutionManager::space() const {
  return IOType::memory(localSolutionLocation_, subgrid_).space();
}

const LocalSolutionManager::LocalGridPartType& LocalSolutionManager::grid_part() const {
  return IOType::memory(localSolutionLocation_, subgrid_).grid_part();
}

void LocalSolutionManager::load() {
  auto& reader = IOType::memory(localSolutionLocation_, subgrid_);
  for (unsigned int i = 0; i < numLocalProblems_; ++i) {
    localSolutions_[i]->clear();
    reader.read(i, localSolutions_[i]);
  }
} // load

void LocalSolutionManager::save() const {
  auto& writer = IOType::memory(localSolutionLocation_, subgrid_);
  for (auto& it : localSolutions_)
    writer.append(it);
} // save

std::size_t LocalSolutionManager::numBoundaryCorrectors() const { return numBoundaryCorrectors_; }

} // namespace MsFEM {
} // namespace Multiscale {
} // namespace Dune {
