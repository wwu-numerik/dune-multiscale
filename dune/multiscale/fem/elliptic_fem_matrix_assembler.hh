// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DiscreteElliptic_HH
#define DiscreteElliptic_HH


#include <memory>
#include <dune/common/fmatrix.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/multiscale/common/traits.hh>

#include <boost/noncopyable.hpp>

namespace Dune {
namespace Multiscale {
namespace FEM {

// used mith multiple type combinations, no real possibilty of splitting

//! \TODO docme
template <class DiscreteFunctionImp, class DiffusionImp>
class DiscreteEllipticOperator
    : public Operator<typename DiscreteFunctionImp::RangeFieldType, typename DiscreteFunctionImp::RangeFieldType,
                      DiscreteFunctionImp, DiscreteFunctionImp>,
      boost::noncopyable {
  typedef DiscreteEllipticOperator<DiscreteFunctionImp, DiffusionImp> This;

  typedef DiscreteFunctionImp DiscreteFunction;

  typedef typename DiscreteFunction::DiscreteFunctionSpaceType DiscreteFunctionSpace;

  typedef typename DiscreteFunctionSpace::GridPartType GridPart;
  typedef typename GridPart::GridType Grid;
  typedef typename DiscreteFunctionSpace::RangeFieldType RangeFieldType;

  typedef typename DiscreteFunctionSpace::DomainType DomainType;
  typedef typename DiscreteFunctionSpace::RangeType RangeType;

  static const int dimension = GridPart::GridType::dimension;
  static const int polynomialOrder = DiscreteFunctionSpace::polynomialOrder;

  typedef typename DiscreteFunctionSpace::BasisFunctionSetType BaseFunctionSet;

public:

  /**
   * \param lower_order_term Operator assumes ownership of it
   **/
  DiscreteEllipticOperator(const DiscreteFunctionSpace& discreteFunctionSpace, const DiffusionImp& diffusion_op,
                           const std::unique_ptr<const CommonTraits::LowerOrderTermType>& lower_order_term = nullptr)
    : discreteFunctionSpace_(discreteFunctionSpace)
    , diffusion_operator_(diffusion_op)
    , lower_order_term_(lower_order_term) {}

private:
  /** dummy implementation of "operator()"
   * \param w = effect of the discrete operator on 'u'
   **/
  virtual void operator()(const DiscreteFunction& u, DiscreteFunction& w) const;

public:
  template <class MatrixType>
  void assemble_matrix(MatrixType& global_matrix) const;

  /** assemble stiffness matrix for the jacobian matrix of the diffusion operator evaluated in the gradient of a certain
   * discrete function (in case of the Newton method, it is the preceeding iterate u_H^{(n-1)} )
   * stiffness matrix with entries
   * \int JA(\nabla disc_func) \nabla phi_i \nabla phi_j
   * (here, JA denotes the jacobian matrix of the diffusion operator A)
   **/
  template <class MatrixType>
  void assemble_jacobian_matrix(DiscreteFunction& disc_func, MatrixType& global_matrix) const;

  // for inhomogeneous boundary condition
  template <class MatrixType>
  void assemble_jacobian_matrix(DiscreteFunction& disc_func, const DiscreteFunction& dirichlet_extension,
                                MatrixType& global_matrix) const;

private:
  const DiscreteFunctionSpace& discreteFunctionSpace_;
  const DiffusionImp& diffusion_operator_;
  const std::unique_ptr<const CommonTraits::LowerOrderTermType>& lower_order_term_;
};

//! \TODO docme
template <class DiscreteFunctionImp, class DiffusionImp>
class SMPDiscreteEllipticOperator : public boost::noncopyable {
  typedef SMPDiscreteEllipticOperator<DiscreteFunctionImp, DiffusionImp> This;

  typedef DiscreteFunctionImp DiscreteFunction;

  typedef typename DiscreteFunction::DiscreteFunctionSpaceType DiscreteFunctionSpace;

  typedef typename DiscreteFunctionSpace::GridPartType GridPart;
  typedef typename GridPart::GridType Grid;
  typedef typename DiscreteFunctionSpace::RangeFieldType RangeFieldType;

  typedef typename DiscreteFunctionSpace::DomainType DomainType;
  typedef typename DiscreteFunctionSpace::RangeType RangeType;

  static const int dimension = GridPart::GridType::dimension;
  static const int polynomialOrder = DiscreteFunctionSpace::polynomialOrder;

  typedef typename DiscreteFunctionSpace::BasisFunctionSetType BaseFunctionSet;

public:

  /**
   * \param lower_order_term Operator assumes ownership of it
   **/
  SMPDiscreteEllipticOperator(const DiscreteFunctionSpace& discreteFunctionSpace, const DiffusionImp& diffusion_op)
    : discreteFunctionSpace_(discreteFunctionSpace)
    , diffusion_operator_(diffusion_op)
    {}

public:
  template <class MatrixType>
  void assemble_matrix(MatrixType& global_matrix) const;

private:
  const DiscreteFunctionSpace& discreteFunctionSpace_;
  const DiffusionImp& diffusion_operator_;
};

} // namespace FEM {
} // namespace Multiscale {
} // namespace Dune {

#include "elliptic_fem_matrix_assembler.cxx"

#endif // #ifndef DiscreteElliptic_HH
