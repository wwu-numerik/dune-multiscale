// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DHUUNE_LOD_TRAITS_HH
#define DHUUNE_LOD_TRAITS_HH

#include <dune/multiscale/common/traits.hh>
#include <dune/subgrid/subgrid.hh>

namespace Dune {
template <int D, class R>
class Subgrid;

namespace Multiscale {
namespace LOD {
class MacroMicroGridSpecifier;
class SubGridList;

// ! type construction for the MSFEM code
struct LODTraits {
  typedef MacroMicroGridSpecifier MacroMicroGridSpecifierType;
  typedef Dune::SubGrid<CommonTraits::GridType::dimension, typename CommonTraits::GridType> SubGridType;
  typedef SubGridList SubGridListType;

  // the following two may change if we intend to use different meshes on coarse and fine level
  typedef typename CommonTraits::GridType::Codim<0>::Entity CoarseEntityType;
  typedef typename CommonTraits::DiscreteFunctionSpaceType::BasisFunctionSetType CoarseBaseFunctionSetType;
};
} // namespace LOD {
} // namespace Multiscale {
} // namespace Dune {}

#endif // DUNE_LOD_TRAITS_HH