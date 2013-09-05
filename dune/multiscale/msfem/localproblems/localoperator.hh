// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef LOCALOPERATOR_HH
#define LOCALOPERATOR_HH

#include <config.h>
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
  typedef MsFEMLocalProblemSolver::SubDiscreteFunctionVectorType SubDiscreteFunctionVectorType;
  typedef CommonTraits::DiffusionType DiffusionOperatorType;
  typedef CommonTraits::NeumannBCType NeumannBoundaryType;

  enum { faceCodim = 1 };

private:
  typedef SubDiscreteFunctionType DiscreteFunction;
  typedef DiffusionOperatorType   DiffusionModel;

  typedef typename DiscreteFunction::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

  typedef typename DiscreteFunctionSpaceType::GridPartType   GridPartType;
  typedef typename DiscreteFunctionSpaceType::GridType       GridType;
  typedef typename DiscreteFunctionSpaceType::RangeFieldType RangeFieldType;

  typedef typename DiscreteFunctionSpaceType::DomainType DomainType;
  typedef typename DiscreteFunctionSpaceType::RangeType  RangeType;
  typedef typename DiscreteFunctionSpaceType::JacobianRangeType
    JacobianRangeType;
  
  static const int dimension = GridPartType::GridType::dimension;
  static const int polynomialOrder = DiscreteFunctionSpaceType::polynomialOrder;

  typedef typename DiscreteFunction::LocalFunctionType LocalFunctionType;

  typedef typename DiscreteFunctionSpaceType::BasisFunctionSetType                   BasisFunctionSetType;
  typedef typename DiscreteFunctionSpaceType::LagrangePointSetType                  LagrangePointSet;
  typedef typename LagrangePointSet::Codim< 1 >::SubEntityIteratorType FaceDofIterator;

  typedef typename DiscreteFunctionSpaceType::IteratorType EntityIteratorType;
  typedef typename EntityIteratorType::Entity                    EntityType;
  typedef typename EntityType::Geometry                    GeometryType;

  typedef typename GridPartType::IntersectionIteratorType IntersectionIterator;
  typedef typename IntersectionIterator::Intersection Intersection;

  typedef Fem::CachingQuadrature< GridPartType, 0 > QuadratureType;

  typedef typename SubGridList::HostDiscreteFunctionType HostDiscreteFunction;
  typedef typename HostDiscreteFunction::DiscreteFunctionSpaceType HostDiscreteFunctionSpaceType;
  typedef typename HostDiscreteFunctionSpaceType::GridPartType HostGridPart;
  typedef typename HostDiscreteFunctionSpaceType::IteratorType HostIterator;
  typedef typename HostIterator::Entity                    HostEntity;
  typedef typename HostEntity::EntityPointer               HostEntityPointer;
  typedef typename HostDiscreteFunction::LocalFunctionType HostLocalFunction;
  typedef typename HostEntity::Geometry                    HostGeometry;
  typedef typename HostGridPart::IntersectionIteratorType  HostIntersectionIterator;
  typedef typename HostDiscreteFunctionSpaceType::LagrangePointSetType HostLagrangePointSet;
  typedef typename HostLagrangePointSet::Codim< faceCodim >::SubEntityIteratorType
    HostGridFaceDofIteratorType;
  typedef MsFEMTraits::CoarseBaseFunctionSetType CoarseBaseFunctionSetType;
  typedef MsFEMTraits::CoarseEntityType CoarseEntityType;
  typedef MsFEMTraits::MacroMicroGridSpecifierType MacroMicroGridSpecifierType;
  
  typedef Fem::CachingQuadrature< HostGridPart, 1 > HostFaceQuadrature;

  typedef typename MsFEMTraits::SubGridListType::SubFaceQuadratureType FaceQuadratureType;
public:

  LocalProblemOperator(const DiscreteFunctionSpaceType& subDiscreteFunctionSpace, const DiffusionModel& diffusion_op);

  //! assemble stiffness matrix for local problems (oversampling strategy 1)
  void assemble_matrix(MsFEMLocalProblemSolver::LocProbFEMMatrixType& global_matrix) const;

  //! assemble stiffness matrix for local problems (oversampling strategy 2 and 3)
  void assemble_matrix(MsFEMLocalProblemSolver::LocProbFEMMatrixType& global_matrix,
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

  void assemble_local_RHS_Dirichlet_corrector(
    const HostDiscreteFunction& dirichlet_extension,
    const SubGridList::CoarseNodeVectorType& coarse_node_vector, /*for constraints*/
    const int& oversampling_strategy,
    // rhs local msfem problem:
    DiscreteFunction& local_problem_RHS) const;

  void assemble_local_RHS_Neumann_corrector(
    const NeumannBoundaryType& neumann_bc,
    const HostDiscreteFunctionSpaceType& host_space,
    const SubGridList::CoarseNodeVectorType& coarse_node_vector, /*for constraints*/
    const int& oversampling_strategy,
    // rhs local msfem problem:
    DiscreteFunction& local_problem_RHS) const;
    
  void assembleAllLocalRHS(const CoarseEntityType& coarseEntity, const MacroMicroGridSpecifierType& specifier,
                  SubDiscreteFunctionVectorType& allLocalRHS) const;

  // assemble various right hand sides (for solving the local saddle point problems with lagrange multpliers)
  void assemble_local_RHS_lg_problems( const HostDiscreteFunction/*CoarseBasisFunctionType*/& coarse_basis_func, double weight,
                                       DiscreteFunction& local_problem_RHS ) const;

  void assemble_local_RHS_lg_problems_all( const std::vector< std::shared_ptr<HostDiscreteFunction > >& coarse_basis_func_list,
                                           std::vector< double >& weights,
                                           std::vector< int >& ids_basis_functions_in_subgrid,
                                           std::vector< std::unique_ptr< DiscreteFunction > >& local_problem_RHS ) const;


  // given a discrete function (representing a right hands side of a local problem,
  // defined on a subgrid) set the boundary dofs to zero
  void set_zero_boundary_condition_RHS(const HostDiscreteFunctionSpaceType& host_space, DiscreteFunction& rhs) const;
  
  void printLocalRHS(const DiscreteFunction& rhs) const;

  double normRHS(const DiscreteFunction& rhs) const;

  // is a given 'point' in the convex hull of corner 0, corner 1 and corner 2 (which forms a codim 0 entity)
  bool point_is_in_element( const DomainType& corner_0, const DomainType& corner_1, const DomainType& corner_2, const DomainType& point) const;
  void projectDirichletValues(HostDiscreteFunction& function) const;
private:
  const DiscreteFunctionSpaceType& subDiscreteFunctionSpace_;
  const DiffusionModel& diffusion_operator_;
};

} //namespace MsFEM {
} //namespace Multiscale {
} //namespace Dune {

#endif // LOCALOPERATOR_HH
