// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef LOCALSOLUTIONMANAGER_HEADERGUARD
#define LOCALSOLUTIONMANAGER_HEADERGUARD

#include <dune/multiscale/common/traits.hh>
#include <dune/multiscale/msfem/msfem_traits.hh>

#include <cstddef>
#include <string>

namespace Dune {
namespace Multiscale {

class MemoryBackend;
class LocalGridList;
/**
 * @brief One LocalSolutionManager instance per coarse cell
 */
class LocalSolutionManager {
private:
  typedef MsFEMTraits::CoarseEntityType CoarseEntityType;
  typedef MsFEMTraits::LocalGridViewType LocalGridViewType;
  typedef MsFEMTraits::LocalGridDiscreteFunctionType LocalGridDiscreteFunctionType;
  typedef MsFEMTraits::LocalSpaceType LocalSpaceType;

public:
  LocalSolutionManager(const CommonTraits::SpaceType& coarse_space, const CoarseEntityType& coarseEntity,
                       const LocalGridList& subgridList);

  MsFEMTraits::LocalSolutionVectorType& getLocalSolutions();

  const LocalSpaceType& space() const;
  const LocalGridViewType& grid_view() const;

  void load();
  void save() const;

  std::size_t numBoundaryCorrectors() const;

private:
  const LocalGridList& subgridList_;
  const MsFEMTraits::LocalGridType& subgrid_;
  LocalGridViewType grid_view_;
  const std::size_t numBoundaryCorrectors_;
  const std::size_t numLocalProblems_;
  MsFEMTraits::LocalSolutionVectorType localSolutions_;
  const std::string localSolutionLocation_;
  MemoryBackend& memory_backend_;
};
}
}

#endif // Header guard
