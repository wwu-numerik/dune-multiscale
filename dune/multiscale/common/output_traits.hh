#ifndef OUTPUT_TRAITS_HH
#define OUTPUT_TRAITS_HH

#include <dune/multiscale/common/traits.hh>
#include <dune/multiscale/problems/elliptic/selector.hh>

namespace Dune {
namespace Multiscale {

struct OutputTraits {
//! --------- typedefs and classes for data output -----------------------------------------
typedef std::tuple< const CommonTraits::DiscreteFunctionType* >      IOTupleType;
typedef Dune::DataOutput< CommonTraits::GridType, IOTupleType > DataOutputType;
typedef Dune::GridFunctionAdapter< Problem::ExactSolution, CommonTraits::GridPartType > DiscreteExactSolutionType;
// just for the discretized exact solution (in case it is available)
typedef std::tuple< const DiscreteExactSolutionType* > ExSolIOTupleType;
// just for the discretized exact solution (in case it is available)
typedef Dune::DataOutput< CommonTraits::GridType, ExSolIOTupleType > ExSolDataOutputType;

};

} // namespace Multiscale
} // namespace Dune

#endif // OUTPUT_TRAITS_HH
