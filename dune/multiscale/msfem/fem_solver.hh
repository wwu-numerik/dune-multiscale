// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef MS_Elliptic_FEM_Solver_HH
#define MS_Elliptic_FEM_Solver_HH


#include <dune/fem/misc/h1norm.hh>
#include <dune/fem/misc/l2error.hh>
#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/solver/cginverseoperator.hh>
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

  typedef CommonTraits::DiscreteFunctionType DiscreteFunctionType;
  typedef DiscreteFunctionType DiscreteFunction;
  typedef typename DiscreteFunction::DiscreteFunctionSpaceType DiscreteFunctionSpace;
  typedef typename DiscreteFunctionSpace::LagrangePointSetType LagrangePointSet;
  typedef typename DiscreteFunctionSpace::GridPartType GridPart;
  typedef typename GridPart::IntersectionIteratorType IntersectionIterator;
  typedef DummyMass<DiscreteFunctionSpace> DummyMassType;

  static const int faceCodim = 1;

private:
  const DiscreteFunctionSpace& discreteFunctionSpace_;


  void solve_linear(const CommonTraits::DiffusionType& diffusion_op,
                    const std::unique_ptr<const CommonTraits::LowerOrderTermType>& lower_order_term,
                    const CommonTraits::FirstSourceType& f, DiscreteFunction& solution, const bool use_smp) const;
  void solve_nonlinear(const CommonTraits::DiffusionType& diffusion_op,
                        const std::unique_ptr<const CommonTraits::LowerOrderTermType>& lower_order_term,
                        const CommonTraits::FirstSourceType& f, DiscreteFunction& solution) const;
public:
  Elliptic_FEM_Solver(const DiscreteFunctionSpace& discreteFunctionSpace);

  //! - ∇ (A(x,∇u)) + b ∇u + c u = f - divG
  //! then:
  //! A --> diffusion operator ('DiffusionOperatorType')
  //! b --> advective part ('AdvectionTermType')
  //! c --> reaction part ('ReactionTermType')
  //! f --> 'first' source term, scalar ('SourceTermType')
  //! G --> 'second' source term, vector valued ('SecondSourceTermType')
  //! this is called from LOD code and therefore mostly copypasta
  void apply_lod(const CommonTraits::DiffusionType& diffusion_op,
             const std::unique_ptr<const CommonTraits::LowerOrderTermType>& lower_order_term,
             const CommonTraits::FirstSourceType& f, const CommonTraits::DiscreteFunctionType& dirichlet_extension,
             const CommonTraits::NeumannBCType& neumann_bc, DiscreteFunction& solution) const;

  void apply(const CommonTraits::DiffusionType& diffusion_op,
             const std::unique_ptr<const CommonTraits::LowerOrderTermType>& lower_order_term,
             const CommonTraits::FirstSourceType& f, DiscreteFunction& solution,
             const bool use_smp = false) const;
};

} // namespace Multiscale
} // namespace Dune

#endif // #ifndef MS_Elliptic_FEM_Solver_HH
