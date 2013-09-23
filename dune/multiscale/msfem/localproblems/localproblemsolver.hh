// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DiscreteEllipticMsFEMLocalProblem_HH
#define DiscreteEllipticMsFEMLocalProblem_HH

#include <config.h>
#include <vector>
#include <memory>

#include <dune/common/fmatrix.hh>
#include <dune/fem/operator/matrix/spmatrix.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/operator/2order/lagrangematrixsetup.hh>
#include <dune/fem/operator/matrix/spmatrix.hh>
#include <dune/fem/gridpart/common/gridpart.hh>
#include <dune/fem/solver/cginverseoperator.hh>
#include <dune/fem/solver/petscsolver.hh>
#include <dune/fem/function/petscdiscretefunction/petscdiscretefunction.hh>

#include <dune/multiscale/common/traits.hh>
#include <dune/multiscale/tools/misc/outputparameter.hh>
#include <dune/multiscale/msfem/localproblems/subgrid-list.hh>
#include <dune/multiscale/msfem/localproblems/weighted-clement-operator.hh>


#include <dune/istl/matrix.hh>

namespace Dune {
namespace Multiscale {
namespace MsFEM {

/** \brief define output parameters for local problems
 *  appends "local_problems" for path
 **/
struct LocalProblemDataOutputParameters
  : public OutputParameters
{
public:
  explicit LocalProblemDataOutputParameters();
};

//! --------------------- the essential local msfem problem solver class ---------------------------
class MsFEMLocalProblemSolver
{
private:
  typedef CommonTraits::DiscreteFunctionType HostDiscreteFunctionType;
  typedef MacroMicroGridSpecifier MacroMicroGridSpecifierType;
  typedef CommonTraits::DiffusionType DiffusionOperatorType;
  typedef CommonTraits::NeumannBCType NeumannBoundaryType;

  //! @todo HostDiscreteFunctionType should be replaced by some kind of coarse function type
  typedef std::vector< std::shared_ptr<HostDiscreteFunctionType> > CoarseBasisFunctionListType;

  //! type of discrete function space
  typedef typename HostDiscreteFunctionType::DiscreteFunctionSpaceType HostDiscreteFunctionSpaceType;
  //! type of (non-discrete )function space
  typedef typename HostDiscreteFunctionSpaceType::FunctionSpaceType FunctionSpaceType;
  typedef typename HostDiscreteFunctionSpaceType::GridPartType HostGridPartType;
  typedef typename HostDiscreteFunctionSpaceType::GridType HostGridType;
  typedef typename HostDiscreteFunctionSpaceType::RangeType RangeType;
  //! type of value of a gradient of a function
  typedef typename HostDiscreteFunctionSpaceType::JacobianRangeType JacobianRangeType;
  typedef typename HostDiscreteFunctionSpaceType::DomainType DomainType;

  typedef typename HostGridType::Traits::LeafIndexSet HostGridLeafIndexSet;
  typedef typename HostGridType::Traits::GlobalIdSet::IdType IdType;
  typedef typename HostDiscreteFunctionSpaceType::IteratorType HostGridEntityIteratorType;
  typedef typename HostGridEntityIteratorType::Entity HostEntityType;
  typedef typename HostEntityType::EntityPointer HostEntityPointerType;
  typedef typename HostGridType::Codim< 0 >::Geometry HostGridEntityGeometry;
  typedef typename HostDiscreteFunctionType::LocalFunctionType HostLocalFunctionType;
  typedef typename HostGridPartType::IntersectionIteratorType HostIntersectionIterator;

  typedef MsFEMTraits::CoarseEntityType CoarseEntityType;
  //! ---------------- typedefs for the SubgridDiscreteFunctionSpace -----------------------
  // ( typedefs for the local grid and the corresponding local ('sub') )discrete space )

  //! type of grid
  typedef typename SubGridList::SubGridType SubGridType;
  //! type of grid part
  typedef typename SubGridList::SubGridPartType SubGridPartType;
  //! type of subgrid discrete function space
  typedef typename SubGridList::SubGridDiscreteFunctionSpaceType SubDiscreteFunctionSpaceType;
  

public:
  //! type of subgrid discrete function
  typedef typename SubGridList::SubGridDiscreteFunctionType SubDiscreteFunctionType;
  typedef std::vector<std::unique_ptr<SubDiscreteFunctionType> > SubDiscreteFunctionVectorType;
private:
  typedef typename SubDiscreteFunctionSpaceType::IteratorType SubgridIteratorType;
  typedef typename SubgridIteratorType::Entity SubgridEntityType;
  typedef typename SubgridEntityType::EntityPointer SubgridEntityPointerType;
  typedef typename SubDiscreteFunctionType::LocalFunctionType SubLocalFunctionType;
  typedef typename SubDiscreteFunctionSpaceType::LagrangePointSetType SGLagrangePointSetType;
  typedef typename SubDiscreteFunctionSpaceType::LagrangePointSetType SubgridLagrangePointSetType;
  typedef typename SubGridType::Codim< 0 >::Geometry SubGridEntityGeometry;
  
  static const int faceCodim = 1;
  typedef typename SubgridLagrangePointSetType::Codim< faceCodim >::SubEntityIteratorType
    SubgridFaceDofIteratorType;

  //! polynomial order of base functions
  enum { polynomialOrder = SubDiscreteFunctionSpaceType::polynomialOrder };

  //! --------------------- istl matrix and vector types -------------------------------------

  //! \TODO diese definitionen machen keinen sinn
  typedef BlockVector< FieldVector< double, 1> > VectorType;
  typedef Matrix< FieldMatrix< double,1,1 > > MatrixType;
  typedef MatrixAdapter< MatrixType, VectorType, VectorType > MatrixOperatorType;
  //typedef SeqGS< MatrixType, VectorType, VectorType > PreconditionerType;
  typedef SeqSOR< MatrixType, VectorType, VectorType > PreconditionerType;
  //typedef BiCGSTABSolver< VectorType > SolverType;
  typedef InverseOperatorResult InverseOperatorResultType;
    
public:
  typedef Dune::Fem::PetscLinearOperator< SubDiscreteFunctionType, SubDiscreteFunctionType > LocProbFEMMatrixType;
private:
  typedef Dune::Fem::PetscInverseOperator< SubDiscreteFunctionType,
                                           LocProbFEMMatrixType >
    InverseLocProbFEMMatrixType;
  
  static std::unique_ptr<InverseLocProbFEMMatrixType> make_inverse_operator(const LocProbFEMMatrixType& problem_matrix);

  typedef WeightedClementOperator WeightedClementOperatorType;
  const HostDiscreteFunctionSpaceType& hostDiscreteFunctionSpace_;
  const DiffusionOperatorType& diffusion_;
  const MacroMicroGridSpecifierType& specifier_;
  SubGridList& subgrid_list_;

  std::vector< std::vector< int > >* ids_relevant_basis_functions_for_subgrid_;
  std::vector< double >* inverse_of_L1_norm_coarse_basis_funcs_;
  const CoarseBasisFunctionListType* coarse_basis_;
  const std::map<int,int>* global_id_to_internal_id_;
  
  const NeumannBoundaryType* neumann_bc_;
  const HostDiscreteFunctionType* dirichlet_extension_;

public:
  /** \brief constructor - with diffusion operator A^{\epsilon}(x)
   * \param subgrid_list cannot be const because Dune::Fem does not provide Gridparts that can be build on a const grid
   **/
  MsFEMLocalProblemSolver(const HostDiscreteFunctionSpaceType& hostDiscreteFunctionSpace,
                          const MacroMicroGridSpecifierType& specifier,
                          SubGridList& subgrid_list,
                          const DiffusionOperatorType& diffusion_operator);

  MsFEMLocalProblemSolver(const HostDiscreteFunctionSpaceType& hostDiscreteFunctionSpace,
                          const MacroMicroGridSpecifierType& specifier,
                          SubGridList& subgrid_list,
                          std::vector< std::vector< int > >& ids_basis_functions_in_subgrid,
                          std::vector< double >& inverse_of_L1_norm_coarse_basis_funcs, // || coarse basis function ||_L1^(-1)
                          const DiffusionOperatorType& diffusion_operator,
                          const CoarseBasisFunctionListType& coarse_basis,
                          const std::map<int,int>& global_id_to_internal_id,
                          const NeumannBoundaryType& neumann_bc,
                          const HostDiscreteFunctionType& dirichlet_extension );

  void solveAllLocalProblems(const CoarseEntityType& coarseCell, SubDiscreteFunctionVectorType& allLocalSolutions) const;

  //! ----------- method: solve the local MsFEM problem ------------------------------------------
  void solvelocalproblem(JacobianRangeType& e,
                         SubDiscreteFunctionType& local_problem_solution,
                         const int coarse_index = -1 ) const;

  // Preprocessing step for the LOD:
  // assemble the two relevant system matrices: one for the corrector problem without contraints
  // and the second of the low dimensional lagrange multiplier (describing the inverse of the schur complement)
  void preprocess_corrector_problems( const int coarse_index,
                                            LocProbFEMMatrixType& locprob_system_matrix,
                                            MatrixType &lm_system_matrix ) const;

  // solve local problem for Local Orthogonal Decomposition Method (LOD) 
  void solve_corrector_problem_lod(JacobianRangeType& e,
                                   LocProbFEMMatrixType& locprob_system_matrix,
                                   MatrixType &lm_system_matrix,
                                   SubDiscreteFunctionType& local_corrector,
                                   const int coarse_index /*= -1*/ ) const;    

  // solve Dirichlet boundary corrector problem for Local Orthogonal Decomposition Method (LOD) 
  void solve_dirichlet_corrector_problem_lod( LocProbFEMMatrixType& locprob_system_matrix,
                                              MatrixType &lm_system_matrix,
                                              SubDiscreteFunctionType& local_corrector,
                                              const int coarse_index /*= -1*/ ) const; 

  // solve Neumann boundary corrector problem for Local Orthogonal Decomposition Method (LOD) 
  void solve_neumann_corrector_problem_lod( LocProbFEMMatrixType& locprob_system_matrix,
                                            MatrixType &lm_system_matrix,
                                            SubDiscreteFunctionType& local_corrector,
                                            const int coarse_index /*= -1*/ ) const; 

  //! create a hostgrid function from a subgridfunction
  void subgrid_to_hostrid_function(const SubDiscreteFunctionType& sub_func,
                                   HostDiscreteFunctionType& host_func) const;

  void output_local_solution(const int coarse_index, const int which,
                             const HostDiscreteFunctionType& host_local_solution) const;

  void output_local_solution(const int coarseIndex, const int which,
                             const SubDiscreteFunctionType& solution) const;


  //! method for solving and saving the solutions of the local msfem problems
  //! for the whole set of macro-entities and for every unit vector e_i
  //! ---- method: solve and save the whole set of local msfem problems -----
  //! Use the host-grid entities of Level 'computational_level' as computational domains for the subgrid computations
  void assemble_all(bool /*silent*/ = true /* state information on subgrids */);

}; // end class

} //namespace MsFEM {
} //namespace Multiscale {
} //namespace Dune {

#endif // #ifndef DiscreteEllipticMsFEMLocalProblem_HH
