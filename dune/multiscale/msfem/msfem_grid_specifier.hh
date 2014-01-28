// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef MSFEM_GRID_SPECIFIER_HH
#define MSFEM_GRID_SPECIFIER_HH

#include <dune/multiscale/common/traits.hh>
#include <cstddef>
#include <vector>

namespace Dune {
namespace Multiscale {
namespace MsFEM {

class MacroMicroGridSpecifier {
  //! \todo DiscreteFunctionSpaceType should be replaced be something like "CoarseDiscreteFunctionSpace" and
  // "FineDiscreteFunctionSpace"
  typedef typename CommonTraits::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
  typedef typename DiscreteFunctionSpaceType::RangeType RangeType;

  static const int faceCodim = 1;

public:
  MacroMicroGridSpecifier();

  //! Get the difference between coarse and fine level
  int getLevelDifference() const;

private:
  // level difference between coarse grid level and fine grid level
  const int coarse_level_fine_level_difference_;
};

} // namespace MsFEM {
} // namespace Multiscale {
} // namespace Dune {

#endif
