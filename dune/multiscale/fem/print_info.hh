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

//! outputs Problem info to output stream
void print_info(const Problem::ProblemContainer& problem, std::ostream& out);
//! write discrete function VTK Output
void write_discrete_function(const DMP::ProblemContainer& problem,
                             CommonTraits::DiscreteFunction_ptr& discrete_solution, const std::string prefix);
void write_discrete_function(CommonTraits::DiscreteFunctionType& discrete_solution, const std::string prefix);

} // namespace Multiscale {
} // namespace Dune {

#endif // PRINT_INFO_HH
