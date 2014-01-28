// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DiscreteEllipticMsFEMLocalProblem_HH
#define DiscreteEllipticMsFEMLocalProblem_HH


#include <dune/common/fmatrix.hh>
#include <dune/common/typetraits.hh>
#include <dune/fem/gridpart/common/gridpart.hh>
#include <dune/fem/operator/2order/lagrangematrixsetup.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/operator/matrix/spmatrix.hh>
#include <dune/fem/operator/matrix/spmatrix.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/solver/cginverseoperator.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/matrix.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solver.hh>
#include <dune/multiscale/common/traits.hh>
#include <dune/multiscale/msfem/localproblems/subgrid-list.hh>
#include <dune/multiscale/msfem/localproblems/weighted-clement-operator.hh>
#include <dune/multiscale/tools/misc/outputparameter.hh>
#include <cstddef>
#include <map>
#include <memory>
#include <vector>

#include "dune/multiscale/common/la_backend.hh"
#include "dune/multiscale/msfem/msfem_grid_specifier.hh"
#include "dune/multiscale/msfem/msfem_traits.hh"
#include "dune/multiscale/problems/base.hh"

namespace Dune {
template <class K, int SIZE> class FieldVector;
}  // namespace Dune

namespace Dune {
namespace Multiscale {
namespace MsFEM {

/** \brief define output parameters for local problems
 *  appends "local_problems" for path
 **/
struct LocalProblemDataOutputParameters : public OutputParameters {
public:
  explicit LocalProblemDataOutputParameters();
};

//! --------------------- the essential local msfem problem solver class ---------------------------
class LocalProblemSolver {
public:
  typedef typename MsFEMTraits::LocalGridDiscreteFunctionType LocalGridDiscreteFunctionType;

private:
  typedef typename MsFEMTraits::LocalGridType LocalGridType;
  typedef typename MsFEMTraits::LocalGridPartType LocalGridPartType;
  typedef typename MsFEMTraits::LocalGridDiscreteFunctionSpaceType LocalGridDiscreteFunctionSpaceType;

  typedef CommonTraits::DiffusionType DiffusionOperatorType;
  typedef CommonTraits::NeumannBCType NeumannBoundaryType;

  //! @todo LocalGridDiscreteFunctionType should be replaced by some kind of coarse function type
  typedef std::vector<std::shared_ptr<LocalGridDiscreteFunctionType>> CoarseBasisFunctionListType;

  typedef typename LocalGridDiscreteFunctionSpaceType::FunctionSpaceType FunctionSpaceType;
  typedef typename LocalGridDiscreteFunctionSpaceType::RangeType RangeType;
  //! type of value of a gradient of a function
  typedef typename LocalGridDiscreteFunctionSpaceType::JacobianRangeType JacobianRangeType;
  typedef typename LocalGridDiscreteFunctionSpaceType::DomainType DomainType;

  typedef typename LocalGridType::Traits::LeafIndexSet LocalGridLeafIndexSet;
  typedef typename LocalGridType::Traits::GlobalIdSet::IdType IdType;
  typedef typename LocalGridDiscreteFunctionSpaceType::IteratorType LocalGridEntityIteratorType;
  typedef typename LocalGridEntityIteratorType::Entity LocalEntityType;
  typedef typename LocalEntityType::EntityPointer LocalEntityPointerType;
  typedef typename LocalGridType::Codim<0>::Geometry LocalGridEntityGeometry;
  typedef typename LocalGridDiscreteFunctionType::LocalFunctionType HostLocalFunctionType;
  typedef typename LocalGridPartType::IntersectionIteratorType HostIntersectionIterator;

  typedef MsFEMTraits::CoarseEntityType CoarseEntityType;
  typedef MsFEMTraits::LocalSolutionVectorType LocalGridDiscreteFunctionVectorType;

  typedef typename LocalGridDiscreteFunctionSpaceType::IteratorType LocalGridIteratorType;
  typedef typename LocalGridIteratorType::Entity LocalGridEntityType;
  typedef typename LocalGridEntityType::EntityPointer LocalGridEntityPointerType;
  typedef typename LocalGridDiscreteFunctionType::LocalFunctionType SubLocalFunctionType;
  typedef typename LocalGridDiscreteFunctionSpaceType::LagrangePointSetType SGLagrangePointSetType;
  typedef typename LocalGridDiscreteFunctionSpaceType::LagrangePointSetType LocalGridLagrangePointSetType;
  typedef typename LocalGridType::Codim<0>::Geometry SubGridEntityGeometry;

  static const int faceCodim = 1;
  typedef typename LocalGridLagrangePointSetType::Codim<faceCodim>::SubEntityIteratorType LocalGridFaceDofIteratorType;

  //! polynomial order of base functions
  static const int polynomialOrder = LocalGridDiscreteFunctionSpaceType::polynomialOrder;

  //! --------------------- istl matrix and vector types -------------------------------------

  //! \TODO diese definitionen machen keinen sinn
  typedef BlockVector<FieldVector<double, 1>> VectorType;
  typedef Matrix<FieldMatrix<double, 1, 1>> MatrixType;
  typedef MatrixAdapter<MatrixType, VectorType, VectorType> MatrixOperatorType;
  // typedef SeqGS< MatrixType, VectorType, VectorType > PreconditionerType;
  typedef SeqSOR<MatrixType, VectorType, VectorType> PreconditionerType;
  // typedef BiCGSTABSolver< VectorType > SolverType;
  typedef InverseOperatorResult InverseOperatorResultType;

public:
  typedef typename BackendChooser<LocalGridDiscreteFunctionSpaceType>::LinearOperatorType LocProbLinearOperatorTypeType;

private:
  typedef typename BackendChooser<LocalGridDiscreteFunctionSpaceType>::InverseOperatorType
  InverseLocProbLinearOperatorTypeType;

  static std::unique_ptr<InverseLocProbLinearOperatorTypeType>
  make_inverse_operator(LocProbLinearOperatorTypeType& problem_matrix);

  typedef WeightedClementOperator WeightedClementOperatorType;
  const DiffusionOperatorType& diffusion_;
  LocalGridList& subgrid_list_;
  const CommonTraits::DiscreteFunctionSpaceType& coarse_space_;

#ifdef ENABLE_LOD_ONLY_CODE
  std::vector<std::vector<std::size_t>>* ids_relevant_basis_functions_for_subgrid_;
  std::vector<double>* inverse_of_L1_norm_coarse_basis_funcs_;
  const CoarseBasisFunctionListType* coarse_basis_;
  const std::map<std::size_t, std::size_t>* global_id_to_internal_id_;

  const NeumannBoundaryType* neumann_bc_;
  const LocalGridDiscreteFunctionType* dirichlet_extension_;
#endif // ENABLE_LOD_ONLY_CODE

public:
  /** \brief constructor - with diffusion operator A^{\epsilon}(x)
   * \param subgrid_list cannot be const because Dune::Fem does not provide Gridparts that can be build on a const grid
   **/
  LocalProblemSolver(const CommonTraits::DiscreteFunctionSpaceType& coarse_space, LocalGridList& subgrid_list,
                          const DiffusionOperatorType& diffusion_operator);

  LocalProblemSolver(const CommonTraits::DiscreteFunctionSpaceType &coarse_space, LocalGridList& subgrid_list, std::vector<std::vector<std::size_t>>& ids_basis_functions_in_subgrid,
      std::vector<double>& inverse_of_L1_norm_coarse_basis_funcs, // || coarse basis function ||_L1^(-1)
      const DiffusionOperatorType& diffusion_operator, const CoarseBasisFunctionListType& coarse_basis,
      const std::map<std::size_t, std::size_t>& global_id_to_internal_id, const NeumannBoundaryType& neumann_bc,
      const LocalGridDiscreteFunctionType& dirichlet_extension);

private:
  void solveAllLocalProblems(const CoarseEntityType& coarseCell,
                             LocalGridDiscreteFunctionVectorType& allLocalSolutions) const;

#ifdef ENABLE_LOD_ONLY_CODE
  void solvelocalproblem(JacobianRangeType& e, LocalGridDiscreteFunctionType& local_problem_solution,
                         const int coarse_index = -1) const;

  // Preprocessing step for the LOD:
  // assemble the two relevant system matrices: one for the corrector problem without contraints
  // and the second of the low dimensional lagrange multiplier (describing the inverse of the schur complement)
  void preprocess_corrector_problems(const int coarse_index, LocProbLinearOperatorTypeType& locprob_system_matrix,
                                     MatrixType& lm_system_matrix) const;

  // solve local problem for Local Orthogonal Decomposition Method (LOD)
  void solve_corrector_problem_lod(JacobianRangeType& e, LocProbLinearOperatorTypeType& locprob_system_matrix,
                                   MatrixType& lm_system_matrix, LocalGridDiscreteFunctionType& local_corrector,
                                   const int coarse_index /*= -1*/) const;

  // solve Dirichlet boundary corrector problem for Local Orthogonal Decomposition Method (LOD)
  void solve_dirichlet_corrector_problem_lod(LocProbLinearOperatorTypeType& locprob_system_matrix,
                                             MatrixType& lm_system_matrix, LocalGridDiscreteFunctionType& local_corrector,
                                             const int coarse_index /*= -1*/) const;

  // solve Neumann boundary corrector problem for Local Orthogonal Decomposition Method (LOD)
  void solve_neumann_corrector_problem_lod(LocProbLinearOperatorTypeType& locprob_system_matrix,
                                           MatrixType& lm_system_matrix, LocalGridDiscreteFunctionType& local_corrector,
                                           const int coarse_index /*= -1*/) const;
#endif // ENABLE_LOD_ONLY_CODE

  void output_local_solution(const int coarse_index, const int which,
                             const LocalGridDiscreteFunctionType& host_local_solution) const;

public:

  /** method for solving and saving the solutions of the local msfem problems
    * for the whole set of macro-entities and for every unit vector e_i
    * ---- method: solve and save the whole set of local msfem problems -----
    * Use the host-grid entities of Level 'computational_level' as computational domains for the subgrid computations
    * **/
  void assembleAndSolveAll(bool /*verbose*/ = false /* state information on subgrids */);

}; // end class

} // namespace MsFEM {
} // namespace Multiscale {
} // namespace Dune {

#endif // #ifndef DiscreteEllipticMsFEMLocalProblem_HH
