// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DiscreteEllipticMsFEMLocalProblem_HH
#define DiscreteEllipticMsFEMLocalProblem_HH

#include <vector>

#include <dune/common/fmatrix.hh>
#include <dune/fem/operator/matrix/spmatrix.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/gridpart/common/gridpart.hh>
#include <dune/fem/operator/2order/lagrangematrixsetup.hh>
#include <dune/fem/solver/cginverseoperator.hh>

#include <dune/multiscale/common/traits.hh>
#include <dune/multiscale/tools/misc/outputparameter.hh>
#include <dune/multiscale/msfem/localproblems/subgrid-list.hh>

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
  typedef typename HostDiscreteFunctionSpaceType::IteratorType HostGridEntityIteratorType;
  typedef typename HostGridEntityIteratorType::Entity HostEntityType;
  typedef typename HostEntityType::EntityPointer HostEntityPointerType;
  typedef typename HostGridType::template Codim< 0 >::Geometry HostGridEntityGeometry;
  typedef typename HostDiscreteFunctionType::LocalFunctionType HostLocalFunctionType;
  typedef typename HostGridPartType::IntersectionIteratorType HostIntersectionIterator;

  //! ---------------- typedefs for the SubgridDiscreteFunctionSpace -----------------------
  // ( typedefs for the local grid and the corresponding local ('sub') )discrete space )

  //! type of grid
  typedef typename SubGridList::SubGridType SubGridType;
  //! type of grid part
  typedef LeafGridPart< SubGridType > SubGridPartType;
  //! type of subgrid discrete function space
  typedef LagrangeDiscreteFunctionSpace< FunctionSpaceType, SubGridPartType, 1 >  // 1=POLORDER
    SubDiscreteFunctionSpaceType;

  //! type of subgrid discrete function
  typedef AdaptiveDiscreteFunction< SubDiscreteFunctionSpaceType > SubDiscreteFunctionType;
  typedef typename SubDiscreteFunctionSpaceType::IteratorType SubgridIteratorType;
  typedef typename SubgridIteratorType::Entity SubgridEntityType;
  typedef typename SubgridEntityType::EntityPointer SubgridEntityPointerType;
  typedef typename SubDiscreteFunctionType::LocalFunctionType SubLocalFunctionType;
  typedef typename SubDiscreteFunctionSpaceType::LagrangePointSetType SGLagrangePointSetType;
  typedef typename SubDiscreteFunctionSpaceType::LagrangePointSetType SubgridLagrangePointSetType;

  enum { faceCodim = 1 };
  typedef typename SubgridLagrangePointSetType::template Codim< faceCodim >::SubEntityIteratorType
    SubgridFaceDofIteratorType;

  //! polynomial order of base functions
  enum { polynomialOrder = SubDiscreteFunctionSpaceType::polynomialOrder };

  struct LocProbMatrixTraits
  {
    typedef SubDiscreteFunctionSpaceType                          RowSpaceType;
    typedef SubDiscreteFunctionSpaceType                          ColumnSpaceType;
    typedef LagrangeMatrixSetup< false >                          StencilType;
    typedef ParallelScalarProduct< SubDiscreteFunctionSpaceType > ParallelScalarProductType;

    template< class M >
    struct Adapter
    {
      typedef LagrangeParallelMatrixAdapter< M > MatrixAdapterType;
    };
  };

  typedef SparseRowMatrixOperator< SubDiscreteFunctionType, SubDiscreteFunctionType,
                                   LocProbMatrixTraits > LocProbFEMMatrix;

  #ifdef SYMMETRIC_DIFFUSION_MATRIX
  typedef Dune::Fem::CGInverseOperator< SubDiscreteFunctionType > InverseLocProbFEMMatrix;
  #else
  // OEMGMRESOp //OEMBICGSQOp // OEMBICGSTABOp
  typedef OEMBICGSTABOp< SubDiscreteFunctionType, LocProbFEMMatrix > InverseLocProbFEMMatrix;
  #endif // ifdef SYMMETRIC_DIFFUSION_MATRIX

  const HostDiscreteFunctionSpaceType& hostDiscreteFunctionSpace_;
  const DiffusionOperatorType& diffusion_;
  const MacroMicroGridSpecifierType& specifier_;
  SubGridList& subgrid_list_;

  const CoarseBasisFunctionListType* coarse_basis_;
  const std::map<int,int>* global_id_to_internal_id_;

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
                          const DiffusionOperatorType& diffusion_operator,
                          const CoarseBasisFunctionListType& coarse_basis,
                          const std::map<int,int>& global_id_to_internal_id );

  //! ----------- method: solve the local MsFEM problem ------------------------------------------
  void solvelocalproblem(JacobianRangeType& e,
                         SubDiscreteFunctionType& local_problem_solution,
                         const int coarse_index = -1 ) const;

  //! create a hostgrid function from a subgridfunction
  void subgrid_to_hostrid_function(const SubDiscreteFunctionType& sub_func,
                                   HostDiscreteFunctionType& host_func);

  void output_local_solution(const int coarse_index, const int which,
                             const HostDiscreteFunctionType& host_local_solution) const;
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
