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

  typedef std::shared_ptr<CommonTraits::GridType> GridPtrType;

public:
  Elliptic_FEM_Solver();

  //! - ∇ (A(x,∇u)) + b ∇u + c u = f - divG
  //! then:
  //! A --> diffusion operator ('DiffusionOperatorType')
  //! b --> advective part ('AdvectionTermType')
  //! c --> reaction part ('ReactionTermType')
  //! f --> 'first' source term, scalar ('SourceTermType')
  CommonTraits::ConstDiscreteFunctionType& solve();

private:
  void apply(CommonTraits::DiscreteFunctionType& solution) const;

  GridPtrType grid_;
  const CommonTraits::SpaceType space_;
  CommonTraits::DiscreteFunctionType solution_;
};

} // namespace Multiscale
} // namespace Dune

#endif // #ifndef MS_Elliptic_FEM_Solver_HH
