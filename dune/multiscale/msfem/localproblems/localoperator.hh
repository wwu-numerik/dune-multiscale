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

#include <dune/multiscale/common/traits.hh>
#include <dune/multiscale/msfem/msfem_traits.hh>
#include <dune/multiscale/msfem/localproblems/localproblemsolver.hh>

namespace Dune {
namespace Multiscale {
namespace MsFEM {

class LocalProblemOperator
{
  typedef MsFEMLocalProblemSolver::SubDiscreteFunctionType SubDiscreteFunctionType;
  typedef CommonTraits::DiffusionType DiffusionOperatorType;

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
  LocalProblemOperator(const DiscreteFunctionSpace& subDiscreteFunctionSpace, const DiffusionModel& diffusion_op);

  //! assemble stiffness matrix for local problems (oversampling strategy 1)
  void assemble_matrix(MsFEMLocalProblemSolver::LocProbFEMMatrix& global_matrix) const;

  //! assemble stiffness matrix for local problems (oversampling strategy 2 and 3)
  void assemble_matrix(MsFEMLocalProblemSolver::LocProbFEMMatrix& global_matrix,
                       const SubGridList::CoarseNodeVectorType& coarse_node_vector /*for constraints*/) const;

  //! assemble the right hand side of a local problem (reconstruction problem on entity)
    //! assemble method for the case of a linear diffusion operator
  //! we compute the following entries for each fine-scale base function phi_h_i:
  //! - \int_{T_0} (A^eps ○ F)(x) ∇ \Phi_H(x_T) · ∇ \phi_h_i(x)
  void assemble_local_RHS(
    // direction 'e'
    const JacobianRangeType& e,
    // rhs local msfem problem:
    DiscreteFunction& local_problem_RHS) const;

  //! assemble method for the case of a linear diffusion operator
  //! in a constraint space, for oversampling strategy 2 and 3
  //! we compute the following entries for each fine-scale base function phi_h_i:
  //! - \int_{T_0} (A^eps ○ F)(x) ∇ \Phi_H(x_T) · ∇ \phi_h_i(x)
  void assemble_local_RHS(
    // direction 'e'
    const JacobianRangeType& e,
    const SubGridList::CoarseNodeVectorType& coarse_node_vector, /*for constraints*/
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

#endif // LOCALOPERATOR_HH
