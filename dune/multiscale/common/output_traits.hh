#ifndef OUTPUT_TRAITS_HH
#define OUTPUT_TRAITS_HH

#include <dune/multiscale/common/traits.hh>

namespace Dune {
namespace Multiscale {

//! typedefs and classes for data output
struct OutputTraits {
  typedef Dune::Fem::GridFunctionAdapter<CommonTraits::ExactSolutionType, CommonTraits::GridPartType>
  DiscreteExactSolutionType;
};

} // namespace Multiscale
} // namespace Dune

#endif // OUTPUT_TRAITS_HH
