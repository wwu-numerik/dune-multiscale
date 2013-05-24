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
#include <dune/multiscale/problems/elliptic/selector.hh>

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

  static const int dimension = GridPart::GridType::dimension;
  
  class OrderedDomainType
  : public DomainType
  {
    public:
  
    OrderedDomainType ( DomainType point ) 
    : DomainType( point )
    {   
    }

    const bool operator< ( OrderedDomainType point_2 ) const
    {

      for (int i = 0; i < dimension; ++i)
      {
        if ( (*this)[i] < point_2[i] ) 
          return true;
        if ( (*this)[i] > point_2[i] ) 
          return false;
      }

      return false;
    }

    const bool operator<= ( OrderedDomainType point_2 ) const
    {
      for (int i = 0; i < dimension; ++i)
      {
        if ( (*this)[i] < point_2[i] ) 
          return true;
        if ( (*this)[i] > point_2[i] ) 
          return false;
      }
      return true;
    }
    
    const bool operator> ( OrderedDomainType point_2 ) const
    {
      for (int i = 0; i < dimension; ++i)
      {
        if ( (*this)[i] > point_2[i] ) 
          return true;
        if ( (*this)[i] < point_2[i] ) 
          return false;
      }
      return false;
    }
    
    const bool operator>= ( OrderedDomainType point_2 ) const
    {
      for (int i = 0; i < dimension; ++i)
      {
        if ( (*this)[i] > point_2[i] ) 
          return true;
        if ( (*this)[i] < point_2[i] ) 
          return false;
      }
      return true;
    }
    
    const bool operator== ( OrderedDomainType point_2 ) const
    {
      for (int i = 0; i < dimension; ++i)
      {
        if ( (*this)[i] > point_2[i] ) 
          return false;
        if ( (*this)[i] < point_2[i] ) 
          return false;
      }
      return true;
    }

  };
  
  // typedef typename HostGrid ::template Codim< 0 > :: template Partition< All_Partition > :: LevelIterator
  // LevelEntityIteratorType;

  typedef typename DiscreteFunctionSpace::IteratorType HostgridIterator;
  typedef HostgridIterator CoarsegridIterator;
  
  typedef typename HostgridIterator::Entity HostEntity;

  typedef typename HostEntity::EntityPointer HostEntityPointer;
  
  typedef typename HostEntity::EntitySeed FineGridEntitySeed;

  // typedef typename HostGrid :: template Codim< 0 > :: template Partition< All_Partition > :: LevelIterator
  // HostGridLevelEntityIterator;

  enum { faceCodim = 1 };

  typedef typename GridPart::IntersectionIteratorType IntersectionIterator;

  // --------------------------- subgrid typedefs ------------------------------------

  typedef MsFEMTraits::SubGridType SubGridType;

  typedef LeafGridPart< SubGridType > SubGridPart;

  // typedef typename SubGridType ::template Codim< 0 > :: template Partition< All_Partition > :: LevelIterator
  // SubGridLevelEntityIteratorType;

  typedef LagrangeDiscreteFunctionSpace< FunctionSpace, SubGridPart, 1 >  // 1=POLORDER
    SubgridDiscreteFunctionSpace;

  typedef AdaptiveDiscreteFunction< SubgridDiscreteFunctionSpace > SubgridDiscreteFunction;

  typedef typename SubgridDiscreteFunctionSpace::IteratorType SubGridIterator;
  
  typedef typename SubGridIterator::Entity SubGridEntity;

  typedef typename SubGridEntity::EntityPointer SubGridEntityPointer;
  
  typedef typename SubGridPart::IntersectionIteratorType  SGIntersectionIterator;


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

   //! for each subgrid, store the vector of basis functions ids that
   //! correspond to interior coarse grid nodes in the subgrid
   // information stored in 'std::vector< std::vector< int > >'
   void assemble_interior_basis_ids( MacroMicroGridSpecifier& specifier,
                                     MsFEMTraits::SubGridListType& subgrid_list,
                                     std::map<int,int>& global_id_to_internal_id,
                                     std::map< OrderedDomainType, int >& coordinates_to_global_coarse_node_id,
                                     std::vector< std::vector< int > >& ids_basis_function_in_subgrid) const;

   void subgrid_to_hostrid_projection( const SubgridDiscreteFunction& sub_func,
                                       DiscreteFunction& host_func) const;

  //! create standard coarse grid basis functions as discrete functions defined on the fine grid
  // ------------------------------------------------------------------------------------
  void add_coarse_basis_contribution( MacroMicroGridSpecifier& specifier,
                                      std::map<int,int>& global_id_to_internal_id,
                                      MsFEMBasisFunctionType& msfem_basis_function_list ) const;


  //! add corrector part to MsFEM basis functions
  void add_corrector_contribution(MacroMicroGridSpecifier &specifier,
                                   std::map<int,int>& global_id_to_internal_id,
                                   MsFEMTraits::SubGridListType& subgrid_list,
                                   MsFEMBasisFunctionType& msfem_basis_function_list ) const;



  template< class DiffusionOperator, class SeedSupportStorage >
  RangeType evaluate_bilinear_form( const DiffusionOperator& diffusion_op,
                                    const DiscreteFunction& func1, const DiscreteFunction& func2,
                                    const SeedSupportStorage& support_of_ms_basis_func_intersection ) const
  {
    RangeType value = 0.0;
#if 1
    int polOrder = 2* DiscreteFunctionSpace::polynomialOrder + 2;
    for (int it_id = 0; it_id < support_of_ms_basis_func_intersection.size(); ++it_id)
    {
      typedef typename HostEntity::template Codim< 0 >::EntityPointer
          HostEntityPointer;

      HostEntityPointer it = discreteFunctionSpace_.grid().entityPointer( support_of_ms_basis_func_intersection[it_id] );
 
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
#endif

//! delete this soon -> just keep it for testing for the moment
#if 0
#if 0
    std::vector< int > ids;
    for (int it_id = 0; it_id < support_of_ms_basis_func_intersection.size(); ++it_id)
    {
      typedef typename HostEntity::template Codim< 0 >::EntityPointer
          HostEntityPointer;

      HostEntityPointer it = discreteFunctionSpace_.grid().entityPointer( support_of_ms_basis_func_intersection[it_id] );
      
      const HostGridLeafIndexSet& coarseGridLeafIndexSet = discreteFunctionSpace_.gridPart().grid().leafIndexSet();
      int id = coarseGridLeafIndexSet.index( *it );
      ids.push_back( id );

    }
#endif  

    int polOrder = 2* DiscreteFunctionSpace::polynomialOrder + 2;
    for (HostgridIterator it = discreteFunctionSpace_.begin(); it != discreteFunctionSpace_.end(); ++it)
    {
      typedef typename HostEntity::template Codim< 0 >::EntityPointer
          HostEntityPointer;

#if 0
      const HostGridLeafIndexSet& coarseGridLeafIndexSet = discreteFunctionSpace_.gridPart().grid().leafIndexSet();
      int id = coarseGridLeafIndexSet.index( *it );
#endif

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

#if 0
         if( std::find( ids.begin(), ids.end(), id ) != ids.end() ) {
	   value += weight * ( diffusive_flux[0] * grad_func_2[0] );}
	   //std::cout << "weight * ( diffusive_flux[0] * grad_func_2[0] ) = " << weight * ( diffusive_flux[0] * grad_func_2[0] ) << std::endl; }
#endif
        value += weight * ( diffusive_flux[0] * grad_func_2[0] );

      }
    }
#endif
    return value;
  }


  // ------------------------------------------------------------------------------------
  template< class DiffusionOperator, class MatrixImp, class SeedSupportStorageList, class RelevantConstellationsList >
  void assemble_matrix( const DiffusionOperator& diffusion_op,
                        MsFEMBasisFunctionType& msfem_basis_function_list_1,
                        MsFEMBasisFunctionType& msfem_basis_function_list_2,
                        SeedSupportStorageList& support_of_ms_basis_func_intersection,
                        RelevantConstellationsList& relevant_constellations,
                        MatrixImp& system_matrix ) const
  {
    for (size_t row = 0; row != system_matrix.N(); ++row)
      for (size_t col = 0; col != system_matrix.M(); ++col)
        system_matrix[row][col] = 0.0;

#ifdef SYMMETRIC_DIFFUSION_MATRIX

   for (unsigned int t = 0; t < relevant_constellations.size(); ++t)
   {
      unsigned int row = get<0>(relevant_constellations[t]);
      unsigned int col = get<1>(relevant_constellations[t]);
      system_matrix[row][col]
        = evaluate_bilinear_form( diffusion_op, *(msfem_basis_function_list_1[row]), *(msfem_basis_function_list_2[col]), support_of_ms_basis_func_intersection[row][col] );
   }

   /* old version without 'relevant_constellations'-vector
   for (size_t row = 0; row != system_matrix.N(); ++row)
    for (size_t col = 0; col <= row; ++col)    
      system_matrix[row][col]
        = evaluate_bilinear_form( diffusion_op, *(msfem_basis_function_list_1[row]), *(msfem_basis_function_list_2[col]), support_of_ms_basis_func_intersection[row][col] );
   */
	
   for (size_t col = 0; col != system_matrix.N(); ++col )
    for (size_t row = 0; row < col; ++row)
      system_matrix[row][col] = system_matrix[col][row];

#else
     
   for (unsigned int t = 0; t < relevant_constellations.size(); ++t)
   {
      unsigned int row = get<0>(relevant_constellations[t]);
      unsigned int col = get<1>(relevant_constellations[t]);
      system_matrix[row][col]
        = evaluate_bilinear_form( diffusion_op, *(msfem_basis_function_list_1[row]), *(msfem_basis_function_list_2[col]), support_of_ms_basis_func_intersection[row][col] );

      if ( row != cole )
      { system_matrix[col][row]
        = evaluate_bilinear_form( diffusion_op, *(msfem_basis_function_list_1[col]), *(msfem_basis_function_list_2[row]), support_of_ms_basis_func_intersection[col][row] ); }
   }

#endif
  }


  // ------------------------------------------------------------------------------------
  template< class SourceTerm, class SeedSupportStorageList, class VectorImp >
  void assemble_rhs( const SourceTerm& f,
                     MsFEMBasisFunctionType& msfem_basis_function_list,
                     SeedSupportStorageList& support_of_ms_basis_func_intersection,
                     VectorImp& rhs ) const
  {

    for (size_t col = 0; col != rhs.N(); ++col)
      rhs[col] = 0.0;

    for (size_t col = 0; col != rhs.N(); ++col)
    {
      
      typedef typename HostEntity::template Codim< 0 >::EntityPointer
          HostEntityPointer;
	  
      const int polOrder = 2* DiscreteFunctionSpace::polynomialOrder + 2;
      for (int it_id = 0; it_id < support_of_ms_basis_func_intersection[col][col].size(); ++it_id)
//      for (const auto& entity : discreteFunctionSpace_)
      {
        HostEntityPointer it = discreteFunctionSpace_.grid().entityPointer( support_of_ms_basis_func_intersection[col][col][it_id] );
        const HostEntity& entity = *it;

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
                            MacroMicroGridSpecifier &specifier,
                            MsFEMTraits::SubGridListType& subgrid_list,
                            DiscreteFunction& coarse_scale_part,
                            DiscreteFunction& fine_scale_part,
                            DiscreteFunction& solution) const;
};

} //namespace MsFEM {
} //namespace Multiscale {
} //namespace Dune {

#endif // #ifndef Elliptic_RIG_MSEM_Solver_HH
