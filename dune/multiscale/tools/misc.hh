// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_MS__TOOLS_MISC_HH
#define DUNE_MS__TOOLS_MISC_HH

#include <dune/gdt/spaces/interface.hh>
#include <dune/multiscale/msfem/msfem_traits.hh>

namespace Dune {
namespace Stuff {
namespace Grid {

inline bool is_simplex_grid(const Multiscale::CommonTraits::DiscreteFunctionSpaceType& space) {
  return space.grid_view()->grid().leafIndexSet().geomTypes(0).size() == 1 &&
         space.grid_view()->grid().leafIndexSet().geomTypes(0)[0].isSimplex();
}

} // namespace Grid
} // namespace Stuff
} // namespace Dune

#endif // DUNE_MS__TOOLS_MISC_HH
