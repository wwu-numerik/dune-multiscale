#ifndef DUNE_MULTISCALE_GRID_CREATION_HH
#define DUNE_MULTISCALE_GRID_CREATION_HH

#include <dune/multiscale/common/traits.hh>

namespace Dune {
namespace Multiscale {

std::pair<std::shared_ptr<CommonTraits::GridType>, std::shared_ptr<CommonTraits::GridType>> make_grids();

} // namespace Multiscale {
} // namespace Dune {

#endif // DUNE_MULTISCALE_GRID_CREATION_HH
