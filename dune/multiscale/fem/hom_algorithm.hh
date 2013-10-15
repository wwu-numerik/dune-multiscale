#ifndef DUNE_MULTISCALE_HOM_ALGORITHM_HH
#define DUNE_MULTISCALE_HOM_ALGORITHM_HH

#include <dune/multiscale/common/traits.hh>

namespace Dune {
namespace Multiscale {
namespace FEM {

//! \TODO docme
void algorithm_hom_fem(typename CommonTraits::GridPointerType& macro_grid_pointer, const std::string filename);

} // namespace FEM {
} // namespace Multiscale {
} // namespace Dune {

#endif // DUNE_MULTISCALE_HOM_ALGORITHM_HH
