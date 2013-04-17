// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef LOCALOPERATOR_HH
#define LOCALOPERATOR_HH

#include <dune/common/fmatrix.hh>
#include <dune/fem/operator/matrix/spmatrix.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/gridpart/common/gridpart.hh>
#include <dune/fem/operator/2order/lagrangematrixsetup.hh>
#include <dune/subgrid/subgrid.hh>
#include <dune/stuff/common/filesystem.hh>
#include <dune/stuff/fem/functions/checks.hh>
#include <dune/multiscale/tools/subgrid_io.hh>
#include <dune/multiscale/tools/disc_func_writer/discretefunctionwriter.hh>
#include <dune/multiscale/tools/misc/outputparameter.hh>
#include <dune/multiscale/msfem/msfem_traits.hh>
#include <dune/multiscale/tools/misc/uzawa.hh>
#include <dune/multiscale/tools/misc/weighted-clement-operator.hh>

namespace Dune {
namespace Multiscale {
namespace MsFEM {

// Imp stands for Implementation
template< class SubDiscreteFunctionType, class DiffusionOperatorType >
class LocalProblemOperator
  : public Operator< typename SubDiscreteFunctionType::RangeFieldType,
                     typename SubDiscreteFunctionType::RangeFieldType,
                     SubDiscreteFunctionType,
                     SubDiscreteFunctionType >
{
  typedef LocalProblemOperator< SubDiscreteFunctionType, DiffusionOperatorType > This;

private:
  typedef SubDiscreteFunctionType DiscreteFunction;
  typedef DiffusionOperatorType   DiffusionModel;

  typedef typename DiscreteFunction::DiscreteFunctionSpaceType DiscreteFunctionSpace;

  typedef typename DiscreteFunctionSpace::GridPartType   GridPart;
  typedef typename DiscreteFunctionSpace::GridType       GridType;
  typedef typename DiscreteFunctionSpace::RangeFieldType RangeFieldType;

  typedef typename DiscreteFunctionSpace::DomainType DomainType;
  typedef typename DiscreteFunctionSpace::RangeType  RangeType;
  typedef typename DiscreteFunctionSpace::JacobianRangeType
    JacobianRangeType;

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
  LocalProblemOperator(const DiscreteFunctionSpace& subDiscreteFunctionSpace, const DiffusionModel& diffusion_op)
    : subDiscreteFunctionSpace_(subDiscreteFunctionSpace)
      , diffusion_operator_(diffusion_op)
  {}

private:
  LocalProblemOperator(const This&) = delete;

public:
  // dummy operator
  virtual void operator()(const DiscreteFunction& u, DiscreteFunction& w) const;

  // assemble stiffness matrix for local problems (oversampling strategy 1)
  template< class MatrixType >
  void assemble_matrix(MatrixType& global_matrix) const;

  // assemble stiffness matrix for local problems (oversampling strategy 2 and 3)
  template< class MatrixType, class CoarseNodeVectorType >
  void assemble_matrix(MatrixType& global_matrix, const CoarseNodeVectorType& coarse_node_vector /*for constraints*/) const;

  // the right hand side assembler methods
  void assemble_local_RHS(
    // direction 'e'
    const JacobianRangeType& e,
    // rhs local msfem problem:
    DiscreteFunction& local_problem_RHS) const;

  // the right hand side assembler methods
  template< class CoarseNodeVectorType >
  void assemble_local_RHS(
    // direction 'e'
    const JacobianRangeType& e,
    const CoarseNodeVectorType& coarse_node_vector, /*for constraints*/
    const int& oversampling_strategy,
    // rhs local msfem problem:
    DiscreteFunction& local_problem_RHS) const;

  void printLocalRHS(const DiscreteFunction& rhs) const;

  double normRHS(const DiscreteFunction& rhs) const;

  // is a given 'point' in the convex hull of corner 0, corner 1 and corner 2 (which forms a codim 0 entity)
  bool point_is_in_element( const DomainType& corner_0, const DomainType& corner_1, const DomainType& corner_2, const DomainType& point) const;

private:
  const DiscreteFunctionSpace& subDiscreteFunctionSpace_;
  const DiffusionModel& diffusion_operator_;
};

} //namespace MsFEM {
} //namespace Multiscale {
} //namespace Dune {

#include "localoperator.cc"

#endif // LOCALOPERATOR_HH
