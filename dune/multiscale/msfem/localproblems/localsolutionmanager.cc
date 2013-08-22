#include <config.h>

// - Dune includes
#include <dune/multiscale/common/traits.hh>
#include <dune/multiscale/msfem/msfem_traits.hh>
#include <dune/multiscale/tools/discretefunctionwriter.hh>
#include <dune/stuff/common/memory.hh>

#include "localsolutionmanager.hh"

namespace Dune {
namespace Multiscale {
namespace MsFEM {
LocalSolutionManager::LocalSolutionManager(const CoarseEntityType& coarseEntity, SubGridListType& subgridList,
                                           const MacroMicroGridSpecifierType& gridSpecifier)
  : subgridList_(subgridList),
    gridSpecifier_(gridSpecifier),
    subGridPart_(subgridList_.getSubGrid(coarseEntity)),
    localDiscreteFunctionSpace_(subGridPart_),
    coarseId_(gridSpecifier_.coarseSpace().gridPart().grid().globalIdSet().id(coarseEntity)),
    loaded_(false),
    numBoundaryCorrectors_((gridSpecifier_.simplexCoarseGrid()) ? 1 : 2),
    numLocalProblems_((gridSpecifier_.simplexCoarseGrid()) ?
                      GridSelector::dimgrid + 1 : gridSpecifier_.coarseSpace().mapper().maxNumDofs() + 2),
    localSolutions_(numLocalProblems_),
    localSolutionLocation_((boost::format("local_problems/_localProblemSolutions_%d")
                            % coarseId_).str())
{
  for (auto& it : localSolutions_)
    it = DSC::make_unique< DiscreteFunctionType >("Local problem Solution", localDiscreteFunctionSpace_);

}

LocalSolutionManager::LocalSolutionVectorType& LocalSolutionManager::getLocalSolutions()
{
  return localSolutions_;
}

const LocalSolutionManager::DiscreteFunctionSpaceType& LocalSolutionManager::getLocalDiscreteFunctionSpace() const
{
  return localDiscreteFunctionSpace_;
}

const LocalSolutionManager::SubGridPartType& LocalSolutionManager::getSubGridPart() const
{
  return subGridPart_;
}

void LocalSolutionManager::loadLocalSolutions()
{
  // reader for the cell problem data file:
  DiscreteFunctionReader reader(localSolutionLocation_);

  for (unsigned int i = 0; i < numLocalProblems_; ++i) {
    localSolutions_[i]->clear();
    reader.read(i, *(localSolutions_[i]));
  }
  loaded_ = true;
  return;
} // loadLocalSolutions

void LocalSolutionManager::saveLocalSolutions() const
{
  // reader for the cell problem data file:
  DiscreteFunctionWriter writer(localSolutionLocation_);

  for (auto& it : localSolutions_)
    writer.append(*it);
} // saveLocalSolutions

bool LocalSolutionManager::solutionsWereLoaded() const
{
  return loaded_;
}

int LocalSolutionManager::numBoundaryCorrectors() const {
  return numBoundaryCorrectors_;
}
}
}
}
