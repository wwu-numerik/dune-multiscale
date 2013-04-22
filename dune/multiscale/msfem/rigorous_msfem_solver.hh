// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef Elliptic_RIG_MSEM_Solver_HH
#define Elliptic_RIG_MSEM_Solver_HH

#include <dune/common/fmatrix.hh>

#include <dune/fem/solver/oemsolver/oemsolver.hh>
#include <dune/fem/solver/inverseoperators.hh>
#include <dune/fem/operator/discreteoperatorimp.hh>
#include <dune/fem/operator/2order/lagrangematrixsetup.hh>
#include <dune/fem/operator/matrix/spmatrix.hh>
#include <dune/fem/space/common/adaptmanager.hh>


#include <dune/istl/matrix.hh>
#include <dune/stuff/fem/functions/checks.hh>
#include <dune/multiscale/common/traits.hh>
#include <dune/multiscale/msfem/msfem_traits.hh>

namespace Dune {
namespace Multiscale {
namespace MsFEM {


class Elliptic_Rigorous_MsFEM_Solver
{
private:
  typedef CommonTraits::DiscreteFunctionType DiscreteFunctionType;
  typedef DiscreteFunctionType DiscreteFunction;

  typedef typename DiscreteFunction::FunctionSpaceType FunctionSpace;

  typedef typename DiscreteFunction::DiscreteFunctionSpaceType DiscreteFunctionSpace;

  typedef typename DiscreteFunction::LocalFunctionType LocalFunction;

  typedef typename DiscreteFunctionSpace::LagrangePointSetType LagrangePointSet;

  typedef typename DiscreteFunctionSpace::GridPartType GridPart;

  typedef CachingQuadrature< GridPart, 0 > CoarseQuadrature;

  typedef typename DiscreteFunctionSpace::GridType HostGrid;

  typedef typename HostGrid::Traits::LeafIndexSet HostGridLeafIndexSet;

  typedef typename HostGrid::Traits::LeafIndexSet CoarseGridLeafIndexSet;

  typedef typename DiscreteFunctionSpace::DomainType        DomainType;
  typedef typename DiscreteFunctionSpace::RangeType         RangeType;
  typedef typename DiscreteFunctionSpace::JacobianRangeType JacobianRangeType;

  // typedef typename HostGrid ::template Codim< 0 > :: template Partition< All_Partition > :: LevelIterator
  // LevelEntityIteratorType;

  typedef typename DiscreteFunctionSpace::IteratorType HostgridIterator;

  typedef typename HostgridIterator::Entity HostEntity;

  typedef typename HostEntity::EntityPointer HostEntityPointer;

  // typedef typename HostGrid :: template Codim< 0 > :: template Partition< All_Partition > :: LevelIterator
  // HostGridLevelEntityIterator;

  enum { faceCodim = 1 };

  typedef typename GridPart::IntersectionIteratorType IntersectionIterator;

  typedef typename LagrangePointSet::template Codim< faceCodim >
    ::SubEntityIteratorType
  FaceDofIterator;

  // --------------------------- subgrid typedefs ------------------------------------

  typedef MsFEMTraits::SubGridType SubGridType;

  typedef LeafGridPart< SubGridType > SubGridPart;

  // typedef typename SubGridType ::template Codim< 0 > :: template Partition< All_Partition > :: LevelIterator
  // SubGridLevelEntityIteratorType;

  typedef LagrangeDiscreteFunctionSpace< FunctionSpace, SubGridPart, 1 >  // 1=POLORDER
    SubgridDiscreteFunctionSpace;

  typedef AdaptiveDiscreteFunction< SubgridDiscreteFunctionSpace > SubgridDiscreteFunction;

  typedef typename SubgridDiscreteFunctionSpace::IteratorType CoarseGridIterator;

  typedef typename CoarseGridIterator::Entity CoarseGridEntity;

  typedef typename CoarseGridEntity::EntityPointer CoarseGridEntityPointer;

  typedef typename SubgridDiscreteFunction::LocalFunctionType CoarseGridLocalFunction;

  typedef typename SubgridDiscreteFunctionSpace::LagrangePointSetType
  CoarseGridLagrangePointSet;

  typedef typename CoarseGridLagrangePointSet::template Codim< faceCodim >
    ::SubEntityIteratorType
  CoarseGridFaceDofIterator;

  //!-----------------------------------------------------------------------------------------

  //! --------------------- istl matrix and vector types -------------------------------------

  typedef BlockVector< FieldVector< double, 1> > VectorType;
  typedef Matrix< FieldMatrix< double,1,1 > > MatrixType;
  typedef MatrixAdapter< MatrixType, VectorType, VectorType > MatrixOperatorType;
  //typedef SeqGS< MatrixType, VectorType, VectorType > PreconditionerType;
  typedef SeqSOR< MatrixType, VectorType, VectorType > PreconditionerType;
  //typedef BiCGSTABSolver< VectorType > SolverType;
  typedef InverseOperatorResult InverseOperatorResultType;

  typedef std::vector< shared_ptr<DiscreteFunction> > MsFEMBasisFunctionType;
  //! ----------------------------------------------------------------------------------------

  const DiscreteFunctionSpace& discreteFunctionSpace_;


public:
  Elliptic_Rigorous_MsFEM_Solver(const DiscreteFunctionSpace& discreteFunctionSpace);

private:

  //! vtk visualization of msfem basis functions
  void vtk_output( MsFEMBasisFunctionType& msfem_basis_function_list,
                   std::string basis_name = "msfem_basis_function" ) const;

  void subgrid_to_hostrid_projection( const SubgridDiscreteFunction& sub_func,
                                      DiscreteFunction& host_func) const;

  //! create standard coarse grid basis functions as discrete functions defined on the fine grid
  // ------------------------------------------------------------------------------------
  void add_coarse_basis_contribution( MacroMicroGridSpecifier< DiscreteFunctionSpace >& specifier,
                                      std::map<int,int>& global_id_to_internal_id,
                                      MsFEMBasisFunctionType& msfem_basis_function_list ) const;


  //! add corrector part to MsFEM basis functions
  void add_corrector_contribution( MacroMicroGridSpecifier< DiscreteFunctionSpace >& specifier,
                                   std::map<int,int>& global_id_to_internal_id,
                                   MsFEMTraits::SubGridListType& subgrid_list,
                                   MsFEMBasisFunctionType& msfem_basis_function_list ) const;



  template< class DiffusionOperator >
  RangeType evaluate_bilinear_form( const DiffusionOperator& diffusion_op, const DiscreteFunction& func1, const DiscreteFunction& func2 ) const
  {
    RangeType value = 0.0;

    int polOrder = 2* DiscreteFunctionSpace::polynomialOrder + 2;
    for (HostgridIterator it = discreteFunctionSpace_.begin(); it != discreteFunctionSpace_.end(); ++it)
    {
      typedef typename HostEntity::template Codim< 0 >::EntityPointer
          HostEntityPointer;

      LocalFunction loc_func_1 = func1.localFunction(*it);
      LocalFunction loc_func_2 = func2.localFunction(*it);

      const auto& geometry = (*it).geometry();

      const CachingQuadrature< GridPart, 0 > quadrature( *it , polOrder);
      const int numQuadraturePoints = quadrature.nop();
      for (int quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint)
      {
        DomainType global_point = geometry.global( quadrature.point(quadraturePoint) );

        //weight
        double weight = geometry.integrationElement( quadrature.point(quadraturePoint) );
        weight *= quadrature.weight(quadraturePoint);

        // gradients of func1 and func2
        JacobianRangeType grad_func_1, grad_func_2;
        loc_func_1.jacobian( quadrature[quadraturePoint], grad_func_1);
        loc_func_2.jacobian( quadrature[quadraturePoint], grad_func_2);

        // A \nabla func1
        JacobianRangeType diffusive_flux(0.0);
        diffusion_op.diffusiveFlux( global_point, grad_func_1, diffusive_flux);

        value += weight * ( diffusive_flux[0] * grad_func_2[0] );

      }
    }
    return value;
  }


  // ------------------------------------------------------------------------------------
  template< class DiffusionOperator, class MatrixImp >
  void assemble_matrix( const DiffusionOperator& diffusion_op,
                        MsFEMBasisFunctionType& msfem_basis_function_list_1,
                        MsFEMBasisFunctionType& msfem_basis_function_list_2,
                        MatrixImp& system_matrix ) const
  {
    for (size_t row = 0; row != system_matrix.N(); ++row)
      for (size_t col = 0; col != system_matrix.M(); ++col)
        system_matrix[row][col] = 0.0;

#ifdef SYMMETRIC_DIFFUSION_MATRIX
   for (size_t row = 0; row != system_matrix.N(); ++row)
    for (size_t col = 0; col <= row; ++col)
      system_matrix[row][col]
        = evaluate_bilinear_form( diffusion_op, *(msfem_basis_function_list_1[row]), *(msfem_basis_function_list_2[col]) );

   for (size_t col = 0; col != system_matrix.N(); ++col )
    for (size_t row = 0; row < col; ++row)
      system_matrix[row][col] = system_matrix[col][row];

#else
   for (size_t row = 0; row != system_matrix.N(); ++row)
     for (size_t col = 0; col != system_matrix.M(); ++col)
       system_matrix[col][row]
         = evaluate_bilinear_form( diffusion_op, *(msfem_basis_function_list_1[row]), *(msfem_basis_function_list_2[col]) );

#endif
  }


  // ------------------------------------------------------------------------------------
  template< class SourceTerm, class VectorImp >
  void assemble_rhs( const SourceTerm& f,
                     MsFEMBasisFunctionType& msfem_basis_function_list,
                     VectorImp& rhs ) const
  {

    for (size_t col = 0; col != rhs.N(); ++col)
      rhs[col] = 0.0;

    for (size_t col = 0; col != rhs.N(); ++col)
    {
      const int polOrder = 2* DiscreteFunctionSpace::polynomialOrder + 2;
      for (const auto& entity : discreteFunctionSpace_)
      {
        const auto& geometry = entity.geometry();

        const auto local_func = msfem_basis_function_list[col]->localFunction(entity);
        const CachingQuadrature< GridPart, 0 > quadrature( entity, polOrder);
        const int numQuadraturePoints = quadrature.nop();
        for (int quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint)
        {
          DomainType global_point = geometry.global( quadrature.point(quadraturePoint) );

          const double weight = geometry.integrationElement( quadrature.point(quadraturePoint) )
                                * quadrature.weight(quadraturePoint);

          // gradients of func1 and func2
          RangeType func_in_x;
          local_func.evaluate( quadrature[quadraturePoint], func_in_x );

          RangeType f_x(0.0);
          f.evaluate( global_point, f_x);

          rhs[col] += weight * ( func_in_x * f_x );
        }
      }
    }
  }

public:
  //! - ∇ (A(x,∇u)) + b ∇u + c u = f - divG
  //! then:
  //! A --> diffusion operator ('DiffusionOperatorType')
  //! b --> advective part ('AdvectionTermType')
  //! c --> reaction part ('ReactionTermType')
  //! f --> 'first' source term, scalar ('SourceTermType')
  //! G --> 'second' source term, vector valued ('SecondSourceTermType')
  //! homogenous Dirchilet boundary condition!:
  void solve_dirichlet_zero(const CommonTraits::DiffusionType& diffusion_op,
                            const CommonTraits::FirstSourceType& f,
                            // number of layers per coarse grid entity T:  U(T) is created by enrichting T with
                            // n(T)-layers.
                            MacroMicroGridSpecifier< DiscreteFunctionSpace >& specifier,
                            MsFEMTraits::SubGridListType& subgrid_list,
                            DiscreteFunction& coarse_scale_part,
                            DiscreteFunction& fine_scale_part,
                            DiscreteFunction& solution) const;
};

} //namespace MsFEM {
} //namespace Multiscale {
} //namespace Dune {

#endif // #ifndef Elliptic_RIG_MSEM_Solver_HH
