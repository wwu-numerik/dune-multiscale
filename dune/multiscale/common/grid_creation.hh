#ifndef DUNE_MULTISCALE_GRID_CREATION_HH
#define DUNE_MULTISCALE_GRID_CREATION_HH

#include <dune/multiscale/common/traits.hh>

namespace Dune {
namespace Multiscale {

//! abstraction for creating coarse and fine grid instances. shared across msfem/fem codes.
std::pair<std::shared_ptr<CommonTraits::GridType>, std::shared_ptr<CommonTraits::GridType>>
make_grids(const bool check_partitioning = true);

} // namespace Multiscale {
} // namespace Dune {

#endif // DUNE_MULTISCALE_GRID_CREATION_HH
