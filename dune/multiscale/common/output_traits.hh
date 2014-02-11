#ifndef OUTPUT_TRAITS_HH
#define OUTPUT_TRAITS_HH

#include <dune/multiscale/common/traits.hh>

namespace Dune {
namespace Multiscale {

//! typedefs and classes for data output
struct OutputTraits {
  typedef std::tuple<const CommonTraits::DiscreteFunctionType*> IOTupleType;
  typedef Dune::Fem::DataWriter<CommonTraits::GridType, IOTupleType> DataOutputType;
  typedef Dune::Fem::GridFunctionAdapter<CommonTraits::ExactSolutionType, CommonTraits::GridPartType>
  DiscreteExactSolutionType;
  typedef std::tuple<const DiscreteExactSolutionType*> ExSolIOTupleType;
  typedef Dune::Fem::DataOutput<CommonTraits::GridType, ExSolIOTupleType> ExSolDataOutputType;
};

} // namespace Multiscale
} // namespace Dune

#endif // OUTPUT_TRAITS_HH
