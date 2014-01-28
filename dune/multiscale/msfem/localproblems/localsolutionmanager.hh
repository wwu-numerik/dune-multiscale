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

namespace Dune {
namespace Multiscale {
namespace MsFEM {

/**
 * @brief One LocalSolutionManager instance per coarse cell
 */
class LocalSolutionManager {
private:
  typedef MsFEMTraits::CoarseEntityType CoarseEntityType;
  typedef MsFEMTraits::LocalGridPartType LocalGridPartType;
  typedef MsFEMTraits::LocalGridDiscreteFunctionType LocalGridDiscreteFunctionType;
  typedef MsFEMTraits::LocalGridDiscreteFunctionSpaceType LocalGridDiscreteFunctionSpaceType;

public:
  LocalSolutionManager(const CommonTraits::DiscreteFunctionSpaceType& coarse_space, const CoarseEntityType& coarseEntity, LocalGridList& subgridList);

  MsFEMTraits::LocalSolutionVectorType& getLocalSolutions();

  const LocalGridDiscreteFunctionSpaceType& space() const;

  const LocalGridPartType& grid_part() const;

  void load();

  void save() const;

  std::size_t numBoundaryCorrectors() const;

private:
  LocalGridList& subgridList_;
  MsFEMTraits::LocalGridType& subgrid_;
  const std::size_t numBoundaryCorrectors_;
  const std::size_t numLocalProblems_;
  MsFEMTraits::LocalSolutionVectorType localSolutions_;
  const std::string localSolutionLocation_;
};
}
}
}

#endif // Header guard
