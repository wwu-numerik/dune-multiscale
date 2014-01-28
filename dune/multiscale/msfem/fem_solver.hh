// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef MS_Elliptic_FEM_Solver_HH
#define MS_Elliptic_FEM_Solver_HH


#include <dune/multiscale/common/traits.hh>
#include <dune/multiscale/fem/fem_traits.hh>
#include <memory>

namespace Dune {
namespace Multiscale {

//! define a dummy mass term:
template <class FunctionSpaceImp>
class DummyMass;

//! \todo docme
class Elliptic_FEM_Solver {

  void solve_linear(const CommonTraits::DiffusionType& diffusion_op,
                    const std::unique_ptr<const CommonTraits::LowerOrderTermType>& lower_order_term,
                    const CommonTraits::FirstSourceType& f, CommonTraits::DiscreteFunctionType& solution, const bool use_smp) const;
  void solve_nonlinear(const CommonTraits::DiffusionType& diffusion_op,
                        const std::unique_ptr<const CommonTraits::LowerOrderTermType>& lower_order_term,
                        const CommonTraits::FirstSourceType& f, CommonTraits::DiscreteFunctionType& solution) const;

  static const int faceCodim = 1;
  const CommonTraits::DiscreteFunctionSpaceType& discreteFunctionSpace_;

public:
  Elliptic_FEM_Solver(const CommonTraits::DiscreteFunctionSpaceType& discreteFunctionSpace);

  //! - ∇ (A(x,∇u)) + b ∇u + c u = f - divG
  //! then:
  //! A --> diffusion operator ('DiffusionOperatorType')
  //! b --> advective part ('AdvectionTermType')
  //! c --> reaction part ('ReactionTermType')
  //! f --> 'first' source term, scalar ('SourceTermType')
  //! G --> 'second' source term, vector valued ('SecondSourceTermType')
  void apply(const CommonTraits::DiffusionType& diffusion_op,
             const std::unique_ptr<const CommonTraits::LowerOrderTermType>& lower_order_term,
             const CommonTraits::FirstSourceType& f, CommonTraits::DiscreteFunctionType& solution,
             const bool use_smp = false) const;
};

} // namespace Multiscale
} // namespace Dune

#endif // #ifndef MS_Elliptic_FEM_Solver_HH
