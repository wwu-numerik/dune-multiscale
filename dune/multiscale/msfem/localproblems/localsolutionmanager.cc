#include <config.h>
#include <boost/format.hpp>
#include <dune/stuff/common/memory.hh>
#include <dune/multiscale/tools/discretefunctionwriter.hh>

#include "localsolutionmanager.hh"

namespace Dune {
namespace Multiscale {
namespace MsFEM {

typedef DiscreteFunctionIO<MsFEMTraits::SubGridDiscreteFunctionType> IOType;

LocalSolutionManager::LocalSolutionManager(const CoarseEntityType& coarseEntity, SubGridListType& subgridList,
                                           const MacroMicroGridSpecifierType& gridSpecifier)
  : subgrid_(subgridList.getSubGrid(coarseEntity))
  , gridSpecifier_(gridSpecifier)
  , coarseId_(gridSpecifier_.coarseSpace().gridPart().grid().globalIdSet().id(coarseEntity))
  , loaded_(false)
  , numBoundaryCorrectors_((gridSpecifier_.simplexCoarseGrid()) ? 1 : 2)
  , numLocalProblems_((gridSpecifier_.simplexCoarseGrid()) ? GridSelector::dimgrid + 1
                                                           : gridSpecifier_.coarseSpace().mapper().maxNumDofs() + 2)
  , localSolutions_(numLocalProblems_)
  , localSolutionLocation_((boost::format("local_problems/_localProblemSolutions_%d") % coarseId_).str())
{
  auto& reader = IOType::memory(localSolutionLocation_, subgrid_);
  for (auto& it : localSolutions_)
    it = make_df_ptr<DiscreteFunctionType>("Local problem Solution", reader.space());
}

const LocalSolutionManager::LocalSolutionVectorType& LocalSolutionManager::getLocalSolutions() const { return localSolutions_; }
LocalSolutionManager::LocalSolutionVectorType& LocalSolutionManager::getLocalSolutions() { return localSolutions_; }

const LocalSolutionManager::DiscreteFunctionSpaceType& LocalSolutionManager::getLocalDiscreteFunctionSpace() const {
  return IOType::memory(localSolutionLocation_, subgrid_).space();
}

const LocalSolutionManager::SubGridPartType& LocalSolutionManager::getSubGridPart() const {
  return IOType::memory(localSolutionLocation_, subgrid_).grid_part();
}

void LocalSolutionManager::loadLocalSolutions() {
  // reader for the cell problem data file:
  auto& reader = IOType::memory(localSolutionLocation_, subgrid_);

  for (unsigned int i = 0; i < numLocalProblems_; ++i) {
    localSolutions_[i]->clear();
    reader.read(i, localSolutions_[i]);
    assert(localSolutions_[i]);
  }
  loaded_ = true;
  return;
} // loadLocalSolutions

void LocalSolutionManager::saveLocalSolutions() const {
  // reader for the cell problem data file:
  auto& writer = IOType::memory(localSolutionLocation_, subgrid_);

  for (const auto& it : localSolutions_)
    writer.append(it);
} // saveLocalSolutions

bool LocalSolutionManager::solutionsWereLoaded() const { return loaded_; }

std::size_t LocalSolutionManager::numBoundaryCorrectors() const { return numBoundaryCorrectors_; }

} // namespace MsFEM
} // namespace Multiscale
} // namespace Dune
