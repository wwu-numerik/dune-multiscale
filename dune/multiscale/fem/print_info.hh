#ifndef PRINT_INFO_HH
#define PRINT_INFO_HH

#include <dune/multiscale/common/traits.hh>
#include <dune/multiscale/problems/base.hh>
#include <dune/multiscale/problems/selector.hh>
#include <dune/stuff/common/ranges.hh>
#include <iosfwd>
#include <string>

namespace Dune {
namespace Multiscale {
namespace FEM {

void print_info(const CommonTraits::ModelProblemDataType& info, std::ostream& out);
void write_discrete_function(CommonTraits::DiscreteFunction_ptr& discrete_solution, const std::string prefix);

} // namespace FEM {
} // namespace Multiscale {
} // namespace Dune {

#endif // PRINT_INFO_HH
