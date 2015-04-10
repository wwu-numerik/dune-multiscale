// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef LOCALOPERATOR_HH
#define LOCALOPERATOR_HH

#include <memory>

#include <dune/multiscale/msfem/msfem_traits.hh>
#include <dune/gdt/assembler/system.hh>
#include <dune/multiscale/problems/base.hh>
#include <dune/multiscale/msfem/diffusion_evaluation.hh>
#include <dune/istl/umfpack.hh>

namespace Dune {
namespace Multiscale {
class LocalProblemOperator {
  typedef typename BackendChooser<MsFEMTraits::LocalSpaceType>::LinearOperatorType LocalLinearOperatorType;
  typedef GDT::Operators::EllipticCG<Problem::LocalDiffusionType, LocalLinearOperatorType, MsFEMTraits::LocalSpaceType>
      EllipticOperatorType;
  typedef GDT::Spaces::Constraints::Dirichlet<typename MsFEMTraits::LocalGridViewType::Intersection,
                                              CommonTraits::RangeFieldType> DirichletConstraintsType;
  typedef DSG::BoundaryInfos::AllDirichlet<MsFEMTraits::LocalGridType::LeafGridView::Intersection> BoundaryInfoType;

  typedef UMFPack<typename LocalLinearOperatorType::BackendType> LocalDirectInverseType;

public:
  LocalProblemOperator(const CommonTraits::SpaceType& coarse_space, const MsFEMTraits::LocalSpaceType& subDiscreteFunctionSpace);

  /** Assemble right hand side vectors for all local problems on one coarse cell.
  *
  * @param[in] coarseEntity The coarse cell.
  * @param[out] allLocalRHS A vector with pointers to the discrete functions for the right hand sides.
  *
  * @note The vector allLocalRHS is assumed to have the correct size and contain pointers to all local rhs
  * functions. The discrete functions in allLocalRHS will be cleared in this function.
  */
  void assemble_all_local_rhs(const MsFEMTraits::CoarseEntityType& coarseEntity, MsFEMTraits::LocalSolutionVectorType& allLocalRHS);

  void apply_inverse(const MsFEMTraits::LocalGridDiscreteFunctionType& current_rhs,
                     MsFEMTraits::LocalGridDiscreteFunctionType& current_solution);

private:
  const MsFEMTraits::LocalSpaceType localSpace_;
  const Problem::LocalDiffusionType local_diffusion_operator_;
  const CommonTraits::SpaceType& coarse_space_;
  LocalLinearOperatorType system_matrix_;
  GDT::SystemAssembler<MsFEMTraits::LocalSpaceType> system_assembler_;
  EllipticOperatorType elliptic_operator_;
  BoundaryInfoType boundaryInfo_;
  DirichletConstraintsType dirichletConstraints_;
  DSG::BoundaryInfos::AllDirichlet<MsFEMTraits::LocalGridType::LeafGridView::Intersection> allLocalDirichletInfo_;
  const bool use_umfpack_;
  std::unique_ptr<LocalDirectInverseType> local_direct_inverse_;
};

} // namespace Multiscale {
} // namespace Dune {

#endif // LOCALOPERATOR_HH
