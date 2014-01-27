// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DHUUNE_LOD_TRAITS_HH
#define DHUUNE_LOD_TRAITS_HH

#include <dune/multiscale/common/traits.hh>
#include <dune/grid/sgrid.hh>

namespace Dune {

namespace Multiscale {
namespace LOD {
class MacroMicroGridSpecifier;
class LocalGridList;

// ! type construction for the MSFEM code
struct LODTraits {
  typedef MacroMicroGridSpecifier MacroMicroGridSpecifier;
  typedef Dune::SGrid<CommonTraits::GridType::dimension, CommonTraits::GridType::dimension> LocalGridType;
  typedef LocalGridList LocalGridList;

  // the following two may change if we intend to use different meshes on coarse and fine level
  typedef typename CommonTraits::GridType::Codim<0>::Entity CoarseEntityType;
  typedef typename CommonTraits::DiscreteFunctionSpaceType::BasisFunctionSetType CoarseBaseFunctionSetType;
};
} // namespace LOD {
} // namespace Multiscale {
} // namespace Dune {}

#endif // DUNE_LOD_TRAITS_HH
