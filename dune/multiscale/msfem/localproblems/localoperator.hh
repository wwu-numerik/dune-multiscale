// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef LOCALOPERATOR_HH
#define LOCALOPERATOR_HH

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
  typedef typename BackendChooser<LocalGridDiscreteFunctionSpaceType>::LinearOperatorType LocalLinearOperatorType;

  static const int dimension = GridPartType::GridType::dimension;

  typedef typename LocalEntityType::EntityPointer LocalEntityPointerType;
  typedef MsFEMTraits::CoarseBaseFunctionSetType CoarseBaseFunctionSetType;
  typedef CommonTraits::DiscreteFunctionSpaceType CoarseSpaceType;
  typedef MsFEMTraits::CoarseEntityType CoarseEntityType;

public:
  LocalProblemOperator(const CoarseSpaceType& coarse_space,
                       const LocalGridDiscreteFunctionSpaceType& subDiscreteFunctionSpace,
                       const DiffusionOperatorType& diffusion_op);

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

  void apply_inverse(const MsFEMTraits::LocalGridDiscreteFunctionType& current_rhs,
                     MsFEMTraits::LocalGridDiscreteFunctionType& current_solution);
private:
  /** Set the dirichlet values to a given discrete function on the sub mesh
  *
  * @param[in, out] function The function in which the values will be set.
  */
  void project_dirichlet_values(CommonTraits::DiscreteFunctionType& function) const;

  //! assemble stiffness matrix for local problems
  void assemble_matrix();

  const LocalGridDiscreteFunctionSpaceType& subDiscreteFunctionSpace_;
  const DiffusionOperatorType& diffusion_operator_;
  const CoarseSpaceType& coarse_space_;
  LocalLinearOperatorType system_matrix_;
  static bool cached_;
  static std::vector<CoarseBaseFunctionSetType::JacobianRangeType> coarseBaseJacs_;
  static std::vector<BasisFunctionSetType::JacobianRangeType> dirichletJacs_;
};

} // namespace MsFEM {
} // namespace Multiscale {
} // namespace Dune {

#endif // LOCALOPERATOR_HH
