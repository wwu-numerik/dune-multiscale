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

#include "dune/multiscale/msfem/localproblems/subgrid-list.hh"
#include "dune/multiscale/msfem/msfem_grid_specifier.hh"
#include "dune/multiscale/problems/base.hh"

namespace Dune {
namespace Multiscale {
namespace MsFEM {

class LocalProblemOperator {
  typedef LocalProblemSolver::LocalGridDiscreteFunctionType LocalGridDiscreteFunctionType;
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
  static const int polynomialOrder = LocalGridDiscreteFunctionSpaceType::polynomialOrder;

  typedef typename LocalEntityType::EntityPointer LocalEntityPointerType;
  typedef MsFEMTraits::CoarseBaseFunctionSetType CoarseBaseFunctionSetType;
  typedef MsFEMTraits::CoarseEntityType CoarseEntityType;

public:
  LocalProblemOperator(const LocalGridDiscreteFunctionSpaceType& subDiscreteFunctionSpace, const DiffusionOperatorType& diffusion_op);

  //! assemble stiffness matrix for local problems (oversampling strategy 1)
  void assemble_matrix(LocalProblemSolver::LocProbLinearOperatorTypeType& global_matrix) const;

  //! assemble stiffness matrix for local problems (oversampling strategy 2 and 3)
  void assemble_matrix(LocalProblemSolver::LocProbLinearOperatorTypeType& global_matrix,
                       const LocalGridList::CoarseNodeVectorType& coarse_node_vector /*for constraints*/) const;

  /** assemble the right hand side of a local problem (reconstruction problem on entity)
   * assemble method for the case of a linear diffusion operator
   * we compute the following entries for each fine-scale base function phi_h_i:
   * - \int_{T_0} (A^eps ○ F)(x) ∇ \Phi_H(x_T) · ∇ \phi_h_i(x)
  **/
  void assemble_local_RHS(
      // direction 'e'
      const JacobianRangeType& e,
      // rhs local msfem problem:
      LocalGridDiscreteFunctionType& local_problem_RHS) const;

  /** assemble method for the case of a linear diffusion operator
   * in a constraint space, for oversampling strategy 2 and 3
   * we compute the following entries for each fine-scale base function phi_h_i:
   * - \int_{T_0} (A^eps ○ F)(x) ∇ \Phi_H(x_T) · ∇ \phi_h_i(x)
   **/
  void assemble_local_RHS(
      // direction 'e'
      const JacobianRangeType& e, const LocalGridList::CoarseNodeVectorType& coarse_node_vector, /*for constraints*/
      const int& oversampling_strategy,
      // rhs local msfem problem:
      LocalGridDiscreteFunctionType& local_problem_RHS) const;

  /** Assemble right hand side vectors for all local problems on one coarse cell.
  *
  * @param[in] coarseEntity The coarse cell.
  * @param[in] specifier A MacroMicroGridSpecifier (needed for access to the coarse base function set).
  * @param[out] allLocalRHS A vector with pointers to the discrete functions for the right hand sides.
  *
  * @note The vector allLocalRHS is assumed to have the correct size and contain pointers to all local rhs
  * functions. The discrete functions in allLocalRHS will be cleared in this function.
  */
  void assembleAllLocalRHS(const CoarseEntityType& coarseEntity, const MacroMicroGridSpecifier& specifier,
                           MsFEMTraits::LocalSolutionVectorType& allLocalRHS) const;

#ifdef ENABLE_LOD_ONLY_CODE
  void assemble_local_RHS_Dirichlet_corrector(
      const LocalGridDiscreteFunctionType& dirichlet_extension,
      const LocalGridList::CoarseNodeVectorType& coarse_node_vector, /*for constraints*/
      const int& oversampling_strategy,
      // rhs local msfem problem:
      LocalGridDiscreteFunctionType& local_problem_RHS) const;

  void
  assemble_local_RHS_Neumann_corrector(const NeumannBoundaryType& neumann_bc,
                                       const LocalGridDiscreteFunctionSpaceType& host_space,
                                       const LocalGridList::CoarseNodeVectorType& coarse_node_vector, /*for constraints*/
                                       const int& oversampling_strategy,
                                       // rhs local msfem problem:
                                       LocalGridDiscreteFunctionType& local_problem_RHS) const;

  // assemble various right hand sides (for solving the local saddle point problems with lagrange multpliers)
  void assemble_local_RHS_lg_problems(const LocalGridDiscreteFunctionType /*CoarseBasisFunctionType*/& coarse_basis_func,
                                      double weight, LocalGridDiscreteFunctionType& local_problem_RHS) const;

  void
  assemble_local_RHS_lg_problems_all(const std::vector<std::shared_ptr<LocalGridDiscreteFunctionType>>& coarse_basis_func_list,
                                     std::vector<double>& weights,
                                     std::vector<std::size_t>& ids_basis_functions_in_subgrid,
                                     std::vector<std::unique_ptr<LocalGridDiscreteFunctionType>>& local_problem_RHS) const;
#endif // ENABLE_LOD_ONLY_CODE

  // given a discrete function (representing a right hands side of a local problem,
  // defined on a subgrid) set the boundary dofs to zero
  void set_zero_boundary_condition_RHS(const LocalGridDiscreteFunctionSpaceType& host_space, LocalGridDiscreteFunctionType& rhs) const;

  double normRHS(const LocalGridDiscreteFunctionType& rhs) const;

  // is a given 'point' in the convex hull of corner 0, corner 1 and corner 2 (which forms a codim 0 entity)
  bool point_is_in_element(const DomainType& corner_0, const DomainType& corner_1, const DomainType& corner_2,
                           const DomainType& point) const;
  void projectDirichletValues(CommonTraits::DiscreteFunctionType &function) const;

private:
  const LocalGridDiscreteFunctionSpaceType& subDiscreteFunctionSpace_;
  const DiffusionOperatorType& diffusion_operator_;
};

} // namespace MsFEM {
} // namespace Multiscale {
} // namespace Dune {

#endif // LOCALOPERATOR_HH
