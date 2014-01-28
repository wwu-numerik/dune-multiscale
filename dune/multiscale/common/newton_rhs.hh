#ifndef DUNE_MULTISCALE_COMMON_NEWTON_RHS_HH
#define DUNE_MULTISCALE_COMMON_NEWTON_RHS_HH


#include <dune/multiscale/common/traits.hh>

namespace Dune {
namespace Multiscale {

//! Assembler for right rand side
//! We assemble the right hand side in a LSE, i.e. f \cdot \Phi_H + G \cdot \nabala \Phi_H
//! we call f the first Source and G the second Source
class NewtonRightHandSide {

public:
/**
 * The rhs-assemble()-methods for non-linear elliptic problems
 * if there is a first source f and a lower order term F:
 * discreteFunction is an output parameter (kind of return value)
 * \param old_u_H from the last iteration step
 **/
static void assemble_for_Newton_method(const CommonTraits::FirstSourceType& f, const CommonTraits::DiffusionType& A,
                                       const CommonTraits::LowerOrderTermType& F,
                                       const  CommonTraits::DiscreteFunctionType& old_u_H,
                                       CommonTraits::DiscreteFunctionType& rhsVector);
};

} // end namespace Multiscale
} // end namespace Dune
#endif // DUNE_MULTISCALE_COMMON_NEWTON_RHS_HH
