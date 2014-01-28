#include <config.h>
#include <assert.h>
#include <dune/common/exceptions.hh>
#include <dune/stuff/common/ranges.hh>
#include <dune/stuff/aliases.hh>
#include <memory>
#include <algorithm>

#include "msfem_grid_specifier.hh"

namespace Dune {
namespace Multiscale {
namespace MsFEM {

MacroMicroGridSpecifier::MacroMicroGridSpecifier()
  : coarse_level_fine_level_difference_(std::numeric_limits<int>::max())
{}

//! Get the difference between coarse and fine level
int MacroMicroGridSpecifier::getLevelDifference() const { return coarse_level_fine_level_difference_; }

} // namespace MsFEM {
} // namespace Multiscale {
} // namespace Dune {
