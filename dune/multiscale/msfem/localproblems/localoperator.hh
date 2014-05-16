// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef LOCALOPERATOR_HH
#define LOCALOPERATOR_HH

#include <dune/multiscale/msfem/msfem_traits.hh>
#include <dune/gdt/assembler/system.hh>

namespace Dune {
namespace Multiscale {
namespace MsFEM {

class LocalProblemOperator {
  typedef MsFEMTraits::LocalGridDiscreteFunctionType LocalGridDiscreteFunctionType;
  typedef CommonTraits::DiffusionType DiffusionOperatorType;
  typedef CommonTraits::NeumannBCType NeumannBoundaryType;

  static const int faceCodim = 1;

  typedef typename MsFEMTraits::LocalGridDiscreteFunctionSpaceType LocalGridDiscreteFunctionSpaceType;
  typedef typename LocalGridDiscreteFunctionSpaceType::BaseFunctionSetType BaseFunctionSetType;
  typedef typename LocalGridDiscreteFunctionSpaceType::EntityType EntityType;
  typedef typename LocalGridDiscreteFunctionSpaceType::EntityType LocalEntityType;
  typedef typename BackendChooser<LocalGridDiscreteFunctionSpaceType>::LinearOperatorType LocalLinearOperatorType;

  typedef typename LocalEntityType::EntityPointer LocalEntityPointerType;
  typedef MsFEMTraits::CoarseBaseFunctionSetType CoarseBaseFunctionSetType;
  typedef CommonTraits::DiscreteFunctionSpaceType CoarseSpaceType;
  typedef MsFEMTraits::CoarseEntityType CoarseEntityType;
  typedef GDT::Operators::EllipticCG< CommonTraits::DiffusionType,
    LocalLinearOperatorType, LocalGridDiscreteFunctionSpaceType > EllipticOperatorType;
  typedef GDT::Constraints::Dirichlet < typename MsFEMTraits::LocalGridViewType::Intersection, CommonTraits::RangeFieldType >
    ConstraintsType;
  typedef DSG::BoundaryInfos::AllDirichlet<MsFEMTraits::LocalGridType::LeafGridView::Intersection> BoundaryInfoType;

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
  //! Compute the number of quadrature points needed for a standard quadrature (needed for
  //! memory reservation for caches
  long getNumQuadPoints(const MsFEMTraits::LocalGridDiscreteFunctionSpaceType& discreteFunctionSpace) const;

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
  GDT::SystemAssembler<LocalGridDiscreteFunctionSpaceType> system_assembler_;
  EllipticOperatorType elliptic_operator_;
  BoundaryInfoType boundaryInfo_;
  ConstraintsType constraints_;
  DMP::ZeroDirichletData dirichlet_;
  static bool cached_;
  static std::vector<CoarseBaseFunctionSetType::JacobianRangeType> coarseBaseJacs_;
  static std::vector<BaseFunctionSetType::JacobianRangeType> dirichletJacs_;
  const bool msfemUsesOversampling_;
};

} // namespace MsFEM {
} // namespace Multiscale {
} // namespace Dune {

#endif // LOCALOPERATOR_HH
