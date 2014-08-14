// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef MS_Elliptic_FEM_Solver_HH
#define MS_Elliptic_FEM_Solver_HH

#include <dune/multiscale/common/traits.hh>
#include <dune/multiscale/problems/base.hh>

namespace Dune {
namespace Multiscale {

//! \todo docme
class Elliptic_FEM_Solver {

  static const int faceCodim = 1;
  const CommonTraits::GdtSpaceType& space_;

public:
  Elliptic_FEM_Solver(const CommonTraits::GdtSpaceType& space);

  //! - ∇ (A(x,∇u)) + b ∇u + c u = f - divG
  //! then:
  //! A --> diffusion operator ('DiffusionOperatorType')
  //! b --> advective part ('AdvectionTermType')
  //! c --> reaction part ('ReactionTermType')
  //! f --> 'first' source term, scalar ('SourceTermType')
  void apply(const Problem::DiffusionBase& diffusion_op,
             const Problem::SourceType& f, CommonTraits::DiscreteFunctionType &solution) const;
};

} // namespace Multiscale
} // namespace Dune

#endif // #ifndef MS_Elliptic_FEM_Solver_HH
