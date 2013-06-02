// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DiscreteElliptic_HH
#define DiscreteElliptic_HH

#include <dune/common/fmatrix.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/operator/common/operator.hh>

#include <boost/noncopyable.hpp>
#include <dune/multiscale/problems/elliptic/selector.hh>

namespace Dune {
namespace Multiscale {
namespace FEM {

// used mith multiple type combinations, no real possibilty of splitting

//! \TODO docme
template< class DiscreteFunctionImp, class DiffusionImp, class LowerOrderTermImp >
class DiscreteEllipticOperator
  : public Operator< typename DiscreteFunctionImp::RangeFieldType, typename DiscreteFunctionImp::RangeFieldType,
                     DiscreteFunctionImp, DiscreteFunctionImp >
    , boost::noncopyable
{
  typedef DiscreteEllipticOperator< DiscreteFunctionImp, DiffusionImp, LowerOrderTermImp > This;

private:
  typedef DiscreteFunctionImp DiscreteFunction;
  typedef DiffusionImp        DiffusionModel;
  typedef LowerOrderTermImp   LowerOrderTerm;

  typedef typename DiscreteFunction::DiscreteFunctionSpaceType DiscreteFunctionSpace;

  typedef typename DiscreteFunctionSpace::GridPartType   GridPart;
  typedef typename GridPart::GridType                    Grid;
  typedef typename DiscreteFunctionSpace::RangeFieldType RangeFieldType;

  typedef typename DiscreteFunctionSpace::DomainType DomainType;
  typedef typename DiscreteFunctionSpace::RangeType  RangeType;

  static const int dimension = GridPart::GridType::dimension;
  static const int polynomialOrder = DiscreteFunctionSpace::polynomialOrder;

  typedef typename DiscreteFunction::LocalFunctionType LocalFunction;

  typedef typename DiscreteFunctionSpace::BaseFunctionSetType                   BaseFunctionSet;
  typedef typename DiscreteFunctionSpace::LagrangePointSetType                  LagrangePointSet;
  typedef typename LagrangePointSet::template Codim< 1 >::SubEntityIteratorType FaceDofIterator;

  typedef typename DiscreteFunctionSpace::IteratorType Iterator;
  typedef typename Iterator::Entity                    Entity;
  typedef typename Entity::Geometry                    Geometry;

  typedef typename GridPart::IntersectionIteratorType IntersectionIterator;
  typedef typename IntersectionIterator::Intersection Intersection;

  typedef CachingQuadrature< GridPart, 0 > Quadrature;

public:

  /**
   * \param lower_order_term Operator assumes ownership of it
   **/
  DiscreteEllipticOperator(const DiscreteFunctionSpace& discreteFunctionSpace,
                           const DiffusionModel& diffusion_op,
                           const std::unique_ptr<const LowerOrderTerm>& lower_order_term = nullptr)
    : discreteFunctionSpace_(discreteFunctionSpace)
    , diffusion_operator_(diffusion_op)
    , lower_order_term_(lower_order_term)
  {}

private:
  /** dummy implementation of "operator()"
   * \param w = effect of the discrete operator on 'u'
   **/
  virtual void operator()(const DiscreteFunction& u, DiscreteFunction& w) const;

public:
  template< class MatrixType >
  void assemble_matrix(MatrixType& global_matrix, bool boundary_treatment = true) const;

  //! Matrix Assembler for local problems on a Subgrid of the Hostgrid:
  template< class MatrixType, class HostDiscreteFunctionSpaceType >
  void assemble_matrix(MatrixType& global_matrix, HostDiscreteFunctionSpaceType& hostSpace,
                       bool boundary_treatment = true) const;

  /** assemble stiffness matrix for the jacobian matrix of the diffusion operator evaluated in the gradient of a certain
   * discrete function (in case of the Newton method, it is the preceeding iterate u_H^{(n-1)} )
   * stiffness matrix with entries
   * \int JA(\nabla disc_func) \nabla phi_i \nabla phi_j
   * (here, JA denotes the jacobian matrix of the diffusion operator A)
   **/
  template< class MatrixType >
  void assemble_jacobian_matrix(DiscreteFunction& disc_func, MatrixType& global_matrix, bool boundary_treatment = true) const;

private:
  const DiscreteFunctionSpace& discreteFunctionSpace_;
  const DiffusionModel& diffusion_operator_;
  const std::unique_ptr<const LowerOrderTermImp>& lower_order_term_;
};

} //namespace FEM {
} //namespace Multiscale {
} //namespace Dune {

#include "elliptic_fem_matrix_assembler.cc"

#endif // #ifndef DiscreteElliptic_HH
