// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef LOCALOPERATOR_HH
#define LOCALOPERATOR_HH

#include <dune/multiscale/msfem/msfem_traits.hh>
#include <dune/gdt/assembler/system.hh>
#include <dune/multiscale/problems/base.hh>
#include <dune/multiscale/msfem/diffusion_evaluation.hh>

namespace Dune {
namespace Multiscale {


class LocalProblemOperator {
  typedef MsFEMTraits::LocalGridDiscreteFunctionType LocalGridDiscreteFunctionType;
  typedef Problem::DiffusionBase DiffusionOperatorType;
  typedef Problem::NeumannBCType NeumannBoundaryType;

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

  typedef GDT::Operators::EllipticCG<Problem::LocalDiffusionType, LocalLinearOperatorType,
                                     LocalGridDiscreteFunctionSpaceType> EllipticOperatorType;
  typedef GDT::Spaces::Constraints::Dirichlet<typename MsFEMTraits::LocalGridViewType::Intersection,
                                      CommonTraits::RangeFieldType> DirichletConstraintsType;
  typedef DSG::BoundaryInfos::AllDirichlet<MsFEMTraits::LocalGridType::LeafGridView::Intersection> BoundaryInfoType;

public:
  LocalProblemOperator(const CoarseSpaceType& coarse_space,
                       const LocalGridDiscreteFunctionSpaceType& subDiscreteFunctionSpace);

  /** Assemble right hand side vectors for all local problems on one coarse cell.
  *
  * @param[in] coarseEntity The coarse cell.
  * @param[out] allLocalRHS A vector with pointers to the discrete functions for the right hand sides.
  *
  * @note The vector allLocalRHS is assumed to have the correct size and contain pointers to all local rhs
  * functions. The discrete functions in allLocalRHS will be cleared in this function.
  */
  void assemble_all_local_rhs(const CoarseEntityType& coarseEntity, MsFEMTraits::LocalSolutionVectorType& allLocalRHS);

  void apply_inverse(const MsFEMTraits::LocalGridDiscreteFunctionType& current_rhs,
                     MsFEMTraits::LocalGridDiscreteFunctionType& current_solution);

private:
  const LocalGridDiscreteFunctionSpaceType& localSpace_;
  const Problem::LocalDiffusionType local_diffusion_operator_;
  const CoarseSpaceType& coarse_space_;
  LocalLinearOperatorType system_matrix_;
  GDT::SystemAssembler<LocalGridDiscreteFunctionSpaceType> system_assembler_;
  EllipticOperatorType elliptic_operator_;
  BoundaryInfoType boundaryInfo_;
  DirichletConstraintsType dirichletConstraints_;
  DSG::BoundaryInfos::AllDirichlet<MsFEMTraits::LocalGridType::LeafGridView::Intersection> allLocalDirichletInfo_;
};


} // namespace Multiscale {
} // namespace Dune {

#endif // LOCALOPERATOR_HH
