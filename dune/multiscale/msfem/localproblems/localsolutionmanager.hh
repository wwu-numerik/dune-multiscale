// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef LOCALSOLUTIONMANAGER_HEADERGUARD
#define LOCALSOLUTIONMANAGER_HEADERGUARD

// - Dune includes
#include <dune/multiscale/common/traits.hh>
#include <dune/multiscale/msfem/localproblems/subgrid-list.hh>
// #include <dune/stuff/fem/functions/checks.hh>


namespace Dune {
namespace Multiscale {
namespace MsFEM {
class LocalSolutionManager
{
private:
  typedef MsFEMTraits::SubGridListType                      SubGridListType;
  typedef typename SubGridListType::SubGridPartType         SubGridPartType;
  typedef typename SubGridListType::SubGridDiscreteFunction DiscreteFunctionType;
  typedef CommonTraits::CoarseEntityType                    CoarseEntityType;
  typedef MsFEMTraits::MacroMicroGridSpecifierType          MacroMicroGridSpecifierType;

public:
  typedef std::vector< std::unique_ptr< DiscreteFunctionType > > LocalSolutionVectorType;
  typedef typename SubGridListType::SubGridDiscreteFunctionSpace DiscreteFunctionSpaceType;

  LocalSolutionManager(const CoarseEntityType& coarseEntity, SubGridListType& subgridList,
                       const MacroMicroGridSpecifierType& gridSpecifier);


  LocalSolutionVectorType& getLocalSolutions();

  const DiscreteFunctionSpaceType& getLocalDiscreteFunctionSpace() const;

  const SubGridPartType& getSubGridPart() const;

  void loadLocalSolutions();

  void saveLocalSolutions() const;

  bool solutionsWereLoaded() const;

private:
  SubGridListType&                   subgridList_;
  const MacroMicroGridSpecifierType& gridSpecifier_;
  SubGridPartType                    subGridPart_;
  DiscreteFunctionSpaceType          localDiscreteFunctionSpace_;
  const int                          coarseId_;
  bool                               loaded_;
  const int                          numLocalProblems_;
  LocalSolutionVectorType            localSolutions_;
  const std::string                  localSolutionLocation_;
};
}
}
}

#endif // Header guard