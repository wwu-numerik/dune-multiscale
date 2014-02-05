// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef LOCALOPERATOR_HH
#define LOCALOPERATOR_HH


#include <dune/common/fmatrix.hh>
#include <dune/fem/gridpart/common/gridpart.hh>
#include <dune/fem/operator/2order/lagrangematrixsetup.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/operator/matrix/spmatrix.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/multiscale/common/traits.hh>
#include <dune/multiscale/msfem/localproblems/localproblemsolver.hh>
#include <dune/multiscale/msfem/msfem_traits.hh>
#include <dune/stuff/common/filesystem.hh>
#include <dune/stuff/fem/functions/checks.hh>
#include <cstddef>
#include <memory>
#include <vector>

#include <dune/multiscale/msfem/localproblems/subgrid-list.hh>
#include <dune/multiscale/problems/base.hh>
#include <dune/multiscale/msfem/msfem_traits.hh>

namespace Dune {
namespace Multiscale {
namespace MsFEM {

class LocalProblemOperator {
  typedef MsFEMTraits::LocalGridDiscreteFunctionType LocalGridDiscreteFunctionType;
  typedef CommonTraits::DiffusionType DiffusionOperatorType;
  typedef CommonTraits::NeumannBCType NeumannBoundaryType;

  static const int faceCodim = 1;

  typedef typename LocalGridDiscreteFunctionType::DiscreteFunctionSpaceType LocalGridDiscreteFunctionSpaceType;
  typedef typename LocalGridDiscreteFunctionSpaceType::GridPartType GridPartType;
  typedef typename LocalGridDiscreteFunctionSpaceType::GridType GridType;
  typedef typename LocalGridDiscreteFunctionSpaceType::RangeFieldType RangeFieldType;
  typedef typename LocalGridDiscreteFunctionSpaceType::DomainType DomainType;
  typedef typename LocalGridDiscreteFunctionSpaceType::RangeType RangeType;
  typedef typename LocalGridDiscreteFunctionSpaceType::JacobianRangeType JacobianRangeType;
  typedef typename LocalGridDiscreteFunctionSpaceType::BasisFunctionSetType BasisFunctionSetType;
  typedef typename LocalGridDiscreteFunctionSpaceType::EntityType EntityType;
  typedef typename LocalGridDiscreteFunctionSpaceType::EntityType LocalEntityType;

  static const int dimension = GridPartType::GridType::dimension;

  typedef typename LocalEntityType::EntityPointer LocalEntityPointerType;
  typedef MsFEMTraits::CoarseBaseFunctionSetType CoarseBaseFunctionSetType;
  typedef CommonTraits::DiscreteFunctionSpaceType CoarseSpaceType;
  typedef MsFEMTraits::CoarseEntityType CoarseEntityType;

public:
  LocalProblemOperator(const CoarseSpaceType &coarse_space, const LocalGridDiscreteFunctionSpaceType& subDiscreteFunctionSpace, const DiffusionOperatorType& diffusion_op);

  //! assemble stiffness matrix for local problems
  void assemble_matrix(LocalProblemSolver::LinearOperatorTypeType& global_matrix) const;

  /** Assemble right hand side vectors for all local problems on one coarse cell.
  *
  * @param[in] coarseEntity The coarse cell.
  * @param[out] allLocalRHS A vector with pointers to the discrete functions for the right hand sides.
  *
  * @note The vector allLocalRHS is assumed to have the correct size and contain pointers to all local rhs
  * functions. The discrete functions in allLocalRHS will be cleared in this function.
  */
  void assemble_all_local_rhs(const CoarseEntityType& coarseEntity,
                              MsFEMTraits::LocalSolutionVectorType& allLocalRHS) const;

  void project_dirichlet_values(CommonTraits::DiscreteFunctionType &function) const;

private:
  const LocalGridDiscreteFunctionSpaceType& subDiscreteFunctionSpace_;
  const DiffusionOperatorType& diffusion_operator_;
  const CoarseSpaceType& coarse_space_;
};

} // namespace MsFEM {
} // namespace Multiscale {
} // namespace Dune {

#endif // LOCALOPERATOR_HH
