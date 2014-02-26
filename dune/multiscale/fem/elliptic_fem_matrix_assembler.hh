// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DiscreteElliptic_HH
#define DiscreteElliptic_HH

#include <memory>
#include <dune/common/fmatrix.hh>

#include<dune/pdelab/localoperator/defaultimp.hh>
#include<dune/pdelab/localoperator/pattern.hh>
#include<dune/pdelab/localoperator/flags.hh>
#include<dune/geometry/referenceelements.hh>
#include<dune/geometry/quadraturerules.hh>

#include <dune/multiscale/common/traits.hh>
#include <dune/multiscale/problems/base.hh>
#include <dune/multiscale/problems/selector.hh>

#include <boost/noncopyable.hpp>

namespace Dune {
namespace Multiscale {
namespace FEM {

//! \TODO docme
class Local_CG_FEM_Operator
    : public Dune::PDELab::NumericalJacobianApplyVolume<Local_CG_FEM_Operator>,
      public Dune::PDELab::NumericalJacobianVolume<Local_CG_FEM_Operator>,
      public Dune::PDELab::FullVolumePattern,
      public Dune::PDELab::LocalOperatorDefaultFlags,
      boost::noncopyable
{
  static const int dimension = CommonTraits::GridType::dimension;

  typedef typename CommonTraits::BasisFunctionSetType BaseFunctionSet;

public:
  // pattern assembly flags
  enum { doPatternVolume = true };

  // residual assembly flags
  enum { doAlphaVolume = true };

  /**
   * \param lower_order_term Operator assumes ownership of it
   **/
  Local_CG_FEM_Operator(const CommonTraits::DiffusionType& diffusion_op, const Dune::Multiscale::CommonTraits::FirstSourceType& source,
                           const std::unique_ptr<const CommonTraits::LowerOrderTermType>& lower_order_term = nullptr)
    : diffusion_operator_(diffusion_op)
    , lower_order_term_(lower_order_term)
    , source_(source)
  {}

  // volume integral depending on test and ansatz functions
  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& /*lfsv*/, R& r) const
  {
    // extract some types
    typedef typename LFSU::Traits::FiniteElementType::
      Traits::LocalBasisType::Traits::DomainFieldType DF;
    typedef typename LFSU::Traits::FiniteElementType::
      Traits::LocalBasisType::Traits::RangeFieldType RF;
    typedef typename LFSU::Traits::FiniteElementType::
      Traits::LocalBasisType::Traits::JacobianType JacobianType;
    typedef typename LFSU::Traits::FiniteElementType::
      Traits::LocalBasisType::Traits::RangeType RangeType;
    typedef typename LFSU::Traits::FiniteElementType::
      Traits::LocalBasisType::Traits::DomainType DomainType;
    typedef typename LFSU::Traits::SizeType size_type;

    const auto& rule = Dune::QuadratureRules<DF,dimension>::rule(eg.geometry().type(),CommonTraits::quadrature_order);
    // loop over quadrature points
    for (auto quad_point : rule ) {
      const auto local_point = quad_point.position();
      const auto global_point = eg.geometry().global(local_point);
      // evaluate basis functions on reference element
      std::vector<RangeType> phi(lfsu.size());
      lfsu.finiteElement().localBasis().evaluateFunction(local_point,phi);

      // evaluate gradient of basis functions on reference element
      std::vector<JacobianType> jacobians(lfsu.size());
      lfsu.finiteElement().localBasis().evaluateJacobian(local_point, jacobians);

      // transform gradients from reference element to real element
      const auto jacobian_inverse = eg.geometry().jacobianInverseTransposed(local_point);
      std::vector<DomainType> gradient_phi(lfsu.size());
      for (size_type i=0; i<lfsu.size(); i++)
        jacobian_inverse.mv(jacobians[i][0],gradient_phi[i]);

      // compute gradient of u
      DomainType gradu(0.0);
      for (size_type i=0; i<lfsu.size(); i++)
        gradu.axpy(x(lfsu,i),gradient_phi[i]);

      // integrate grad u * A * grad phi_i - (F+f) phi_i
      Multiscale::Problem::JacobianRangeType gradient_phi_fem;
      Multiscale::Problem::JacobianRangeType diffusion_in_gradient_phi;
      const auto factor = quad_point.weight()*eg.geometry().integrationElement(local_point);
      for (size_type i=0; i<lfsu.size(); i++) {
        gradient_phi_fem[0] = gradient_phi[i];
        diffusion_operator_.diffusiveFlux(global_point, gradient_phi_fem, diffusion_in_gradient_phi);
        RangeType F_i(0.0);
        RangeType f_i(0.0);
        if (lower_order_term_)
          lower_order_term_->evaluate(global_point, phi[i], gradient_phi_fem, F_i);
        source_.evaluate(global_point, f_i);
        F_i += f_i;
        r.accumulate(lfsu, i, ( gradu*diffusion_in_gradient_phi[0] - F_i*phi[i] )*factor);
      }
    }
  }
#if 0
  /** assemble stiffness matrix for the jacobian matrix of the diffusion operator evaluated in the gradient of a certain
   * discrete function (in case of the Newton method, it is the preceeding iterate u_H^{(n-1)} )
   * stiffness matrix with entries
   * \int JA(\nabla disc_func) \nabla phi_i \nabla phi_j
   * (here, JA denotes the jacobian matrix of the diffusion operator A)
   **/
  void assemble_jacobian_matrix(DiscreteFunction& disc_func, MatrixType& global_matrix) const;

  //! for inhomogeneous boundary condition
  void assemble_jacobian_matrix(DiscreteFunction& disc_func, const DiscreteFunction& dirichlet_extension,
                                MatrixType& global_matrix) const;
#endif //0

private:
  const CommonTraits::DiffusionType& diffusion_operator_;
  const std::unique_ptr<const CommonTraits::LowerOrderTermType>& lower_order_term_;
  const Dune::Multiscale::CommonTraits::FirstSourceType& source_;
};


} // namespace FEM {
} // namespace Multiscale {
} // namespace Dune {

#endif // #ifndef DiscreteElliptic_HH
