// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_MS__TOOLS_MISC_HH
#define DUNE_MS__TOOLS_MISC_HH

#include <dune/fem/space/common/discretefunctionspace.hh>
#include <dune/multiscale/msfem/msfem_traits.hh>

namespace Dune {
namespace Stuff {
namespace Grid {

template < class TraitsImp >
bool is_simplex_grid(const Dune::Fem::DiscreteFunctionSpaceInterface<TraitsImp>& space) {
  return space.grid_part().grid().leafIndexSet().geomTypes(0).size() == 1 &&
         space.grid_part().grid().leafIndexSet().geomTypes(0)[0].isSimplex();
}

} // namespace Grid
} // namespace Stuff
} // namespace Dune

#endif // DUNE_MS__TOOLS_MISC_HH
