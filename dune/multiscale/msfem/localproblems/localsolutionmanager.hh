// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef LOCALSOLUTIONMANAGER_HEADERGUARD
#define LOCALSOLUTIONMANAGER_HEADERGUARD

// - Dune includes
#include <dune/multiscale/common/traits.hh>
#include <dune/multiscale/msfem/localproblems/subgrid-list.hh>
#include <cstddef>
#include <memory>
#include <string>
#include <vector>

#include "dune/multiscale/msfem/msfem_grid_specifier.hh"
#include "dune/multiscale/msfem/msfem_traits.hh"
// #include <dune/stuff/fem/functions/checks.hh>

namespace Dune {
namespace Multiscale {
namespace MsFEM {

/**
 * @brief One LocalSolutionManager instance per
 */
class LocalSolutionManager {
private:
  typedef MsFEMTraits::LocalGridListType LocalGridListType;
  typedef typename CommonTraits::GridType::Traits::GlobalIdSet::IdType IdType;
  typedef MsFEMTraits::CoarseEntityType CoarseEntityType;


  typedef MsFEMTraits::LocalGridPartType LocalGridPartType;
  typedef MsFEMTraits::LocalGridDiscreteFunctionType LocalGridDiscreteFunctionType;
  typedef MsFEMTraits::LocalGridDiscreteFunctionSpaceType LocalGridDiscreteFunctionSpaceType;
  typedef MsFEMTraits::MacroMicroGridSpecifierType MacroMicroGridSpecifierType;
public:
  typedef std::vector<std::unique_ptr<LocalGridDiscreteFunctionType>> LocalSolutionVectorType;

  LocalSolutionManager(const CoarseEntityType& coarseEntity, LocalGridListType& subgridList,
                       const MacroMicroGridSpecifierType& gridSpecifier);

  LocalSolutionVectorType& getLocalSolutions();

  const LocalGridDiscreteFunctionSpaceType& getLocalDiscreteFunctionSpace() const;

  const LocalGridPartType& getSubGridPart() const;

  void loadLocalSolutions();

  void saveLocalSolutions() const;

  bool solutionsWereLoaded() const;

  std::size_t numBoundaryCorrectors() const;

private:
  LocalGridListType& subgridList_;
  const MacroMicroGridSpecifierType& gridSpecifier_;
  LocalGridPartType subGridPart_;
  LocalGridDiscreteFunctionSpaceType localDiscreteFunctionSpace_;
  const IdType coarseId_;
  bool loaded_;
  const std::size_t numBoundaryCorrectors_;
  const std::size_t numLocalProblems_;
  LocalSolutionVectorType localSolutions_;
  const std::string localSolutionLocation_;
};
}
}
}

#endif // Header guard
