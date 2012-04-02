#ifndef DiscreteEllipticMSFEMOperator_HH
#define DiscreteEllipticMSFEMOperator_HH

#include <dune/common/fmatrix.hh>

#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/multiscale/tools/solver/MsFEM/msfem_localproblems/subgrid-list.hh>
#include <dune/multiscale/tools/solver/MsFEM/msfem_localproblems/localproblemsolver.hh>

#include <dune/fem/operator/2order/lagrangematrixsetup.hh>

namespace Dune
{  
  // Imp stands for Implementation
  template< class CoarseDiscreteFunctionImp, class MacroMicroGridSpecifierImp, class FineDiscreteFunctionImp, class DiffusionImp >
  class DiscreteEllipticMsFEMOperator
  : public Operator< typename CoarseDiscreteFunctionImp::RangeFieldType,
                     typename CoarseDiscreteFunctionImp::RangeFieldType,
		      CoarseDiscreteFunctionImp, CoarseDiscreteFunctionImp >
  {
    typedef DiscreteEllipticMsFEMOperator< CoarseDiscreteFunctionImp, MacroMicroGridSpecifierImp, FineDiscreteFunctionImp, DiffusionImp > This;

  public:
    
    typedef CoarseDiscreteFunctionImp CoarseDiscreteFunction;
    typedef FineDiscreteFunctionImp FineDiscreteFunction;
    typedef MacroMicroGridSpecifierImp MacroMicroGridSpecifierType;
   
    typedef DiffusionImp DiffusionModel;
        
    typedef typename CoarseDiscreteFunction::DiscreteFunctionSpaceType CoarseDiscreteFunctionSpace;    
    typedef typename FineDiscreteFunction::DiscreteFunctionSpaceType FineDiscreteFunctionSpace;
    
    typedef typename FineDiscreteFunctionSpace :: FunctionSpaceType FunctionSpace;

    typedef typename FineDiscreteFunctionSpace::GridPartType FineGridPart;
    typedef typename FineDiscreteFunctionSpace::GridType FineGrid;

    typedef typename FineDiscreteFunctionSpace::RangeFieldType RangeFieldType;
    typedef typename FineDiscreteFunctionSpace::DomainType DomainType;
    typedef typename FineDiscreteFunctionSpace::RangeType RangeType;
    typedef typename FineDiscreteFunctionSpace::JacobianRangeType
      JacobianRangeType;
      

    typedef SubGrid< WORLDDIM , FineGrid > SubGridType;
    typedef SubGridList< FineDiscreteFunction, SubGridType, MacroMicroGridSpecifierType > SubGridListType;
    
    typedef MsFEMLocalProblemSolver< FineDiscreteFunction, SubGridListType, DiffusionModel > MsFEMLocalProblemSolverType;
    

  protected:
    static const int dimension = FineGridPart::GridType::dimension;
    static const int polynomialOrder = FineDiscreteFunctionSpace::polynomialOrder;

    typedef typename FineDiscreteFunction::LocalFunctionType FineLocalFunction;

    typedef typename FineDiscreteFunctionSpace::BaseFunctionSetType FineBaseFunctionSet;
    typedef typename FineDiscreteFunctionSpace::LagrangePointSetType FineLagrangePointSet;
    typedef typename FineLagrangePointSet::template Codim< 1 >::SubEntityIteratorType FineFaceDofIterator;

    typedef typename FineDiscreteFunctionSpace::IteratorType FineIterator;
    typedef typename FineIterator::Entity FineEntity;
    typedef typename FineEntity::EntityPointer FineEntityPointer; 
    typedef typename FineEntity::Geometry FineGeometry;

    typedef typename FineGridPart::IntersectionIteratorType FineIntersectionIterator;
    typedef typename FineIntersectionIterator::Intersection FineIntersection;

    typedef typename FineGrid :: template Codim< 0 > :: template  Partition< All_Partition > :: LevelIterator FineLevelEntityIterator;

    typedef CachingQuadrature< FineGridPart, 0 > FineQuadrature;


  public:

    typedef typename CoarseDiscreteFunctionSpace::GridPartType CoarseGridPart;
    typedef typename CoarseDiscreteFunctionSpace::GridType CoarseGrid;

  protected:

    typedef typename CoarseDiscreteFunction::LocalFunctionType CoarseLocalFunction;

    typedef typename CoarseDiscreteFunctionSpace::BaseFunctionSetType CoarseBaseFunctionSet;
    typedef typename CoarseDiscreteFunctionSpace::LagrangePointSetType CoarseLagrangePointSet;
    typedef typename CoarseLagrangePointSet::template Codim< 1 >::SubEntityIteratorType CoarseFaceDofIterator;

    typedef typename CoarseDiscreteFunctionSpace::IteratorType CoarseIterator;
    typedef typename CoarseIterator::Entity CoarseEntity;
    typedef typename CoarseEntity::Geometry CoarseGeometry;

    typedef typename CoarseGridPart::IntersectionIteratorType CoarseIntersectionIterator;
    typedef typename CoarseIntersectionIterator::Intersection CoarseIntersection;

    typedef typename CoarseGrid :: template Codim< 0 > :: template  Partition< All_Partition > :: LevelIterator CoarseLevelEntityIterator;

    typedef CachingQuadrature< CoarseGridPart, 0 > CoarseQuadrature;
    
  public:
    
    //! ---------------- typedefs for the SubgridDiscreteFunctionSpace -----------------------
    //  ( typedefs for the local grid and the corresponding local ('sub') )discrete space ) 

    //! type of grid part
    typedef LeafGridPart< SubGridType > SubGridPart; 

    //! type of subgrid discrete function space
    typedef LagrangeDiscreteFunctionSpace < FunctionSpace, SubGridPart, 1 > //1=POLORDER
          LocalDiscreteFunctionSpace;

    //! type of subgrid discrete function
    typedef AdaptiveDiscreteFunction < LocalDiscreteFunctionSpace > LocalDiscreteFunction;

    typedef typename LocalDiscreteFunctionSpace :: IteratorType LocalGridIterator;
    
    typedef typename LocalGridIterator :: Entity LocalGridEntity;
    
    typedef typename LocalGridEntity :: EntityPointer LocalGridEntityPointer;
    
    typedef typename LocalDiscreteFunction :: LocalFunctionType LocalGridLocalFunction;
    
    typedef typename LocalDiscreteFunctionSpace :: LagrangePointSetType LGLagrangePointSet;
    
    typedef typename LocalDiscreteFunctionSpace :: BaseFunctionSetType LocalGridBaseFunctionSet;
    
    typedef typename LocalGridEntity::Geometry LocalGridGeometry;
    
    typedef CachingQuadrature< SubGridPart, 0 > LocalGridQuadrature;

//!-----------------------------------------------------------------------------------------

  public:

    DiscreteEllipticMsFEMOperator( const CoarseDiscreteFunctionSpace &coarseDiscreteFunctionSpace,
                                   const FineDiscreteFunctionSpace &fineDiscreteFunctionSpace,
                                   // number of layers per coarse grid entity T:  U(T) is created by enrichting T with n(T)-layers:
                                   SubGridListType& subgrid_list,
                                   const DiffusionModel &diffusion_op,
                                   std :: ofstream& data_file,
                                   std :: string path = ""  )
    : coarseDiscreteFunctionSpace_( coarseDiscreteFunctionSpace ),
      fineDiscreteFunctionSpace_( fineDiscreteFunctionSpace ),
      subgrid_list_( subgrid_list ),
      diffusion_operator_( diffusion_op ),
      data_file_( &data_file ),
      path_( path )
    {

      bool silence = false;

      const int coarse_level = coarseDiscreteFunctionSpace_.gridPart().grid().maxLevel();

      std :: string local_path =  path_ + "/local_problems/";

      MsFEMLocalProblemSolverType loc_prob_solver( fineDiscreteFunctionSpace_, subgrid_list_, diffusion_operator_, data_file, local_path );
      loc_prob_solver.assemble_all( coarse_level, silence );
    }


  private:
    DiscreteEllipticMsFEMOperator ( const This & );

  public:

    // dummy operator
    virtual void
    operator() ( const CoarseDiscreteFunction &u, CoarseDiscreteFunction &w ) const;


    template< class MatrixType >
    void assemble_matrix ( MatrixType &global_matrix ) const;

  private:

    // create a hostgrid function from a subgridfunction
    // Note: the maximum gride levels for both underlying grids must be the same
    void subgrid_to_hostrid_function( const LocalDiscreteFunction &sub_func, FineDiscreteFunction &host_func );

  private:

    const FineDiscreteFunctionSpace &fineDiscreteFunctionSpace_;

    const CoarseDiscreteFunctionSpace &coarseDiscreteFunctionSpace_;
     
    const DiffusionModel &diffusion_operator_;

    SubGridListType& subgrid_list_;
    
    // data file for saving information
    std :: ofstream *data_file_;

    // path where to save the data output
    std :: string path_;

  };

  // create a hostgrid function from a subgridfunction
  template< class CoarseDiscreteFunctionImp, class MacroMicroGridSpecifierImp, class FineDiscreteFunctionImp, class DiffusionImp >
  void DiscreteEllipticMsFEMOperator< CoarseDiscreteFunctionImp, 
                                      MacroMicroGridSpecifierImp,
                                      FineDiscreteFunctionImp, 
                                      DiffusionImp > :: 
  subgrid_to_hostrid_function( const LocalDiscreteFunction &sub_func, FineDiscreteFunction &host_func )
    {

       if ( sub_func.space().gridPart().grid().maxLevel() != host_func.space().gridPart().grid().maxLevel() )
         { std :: cout << "Error in method 'subgrid_to_hostrid_function': MaxLevel of SubGrid not identical to MaxLevel of FineGrid." << std :: endl; }

       host_func.clear();

       const LocalDiscreteFunctionSpace &subDiscreteFunctionSpace = sub_func.space();
       const SubGridType &subGrid = subDiscreteFunctionSpace.grid();

       LocalGridIterator sub_endit = subDiscreteFunctionSpace.end();
       for( LocalGridIterator sub_it = subDiscreteFunctionSpace.begin(); sub_it != sub_endit; ++sub_it )
          {

             const LocalGridEntity &sub_entity = *sub_it;

             FineEntityPointer host_entity_pointer = subGrid.template getHostEntity<0>( *sub_it );
             const FineEntity& host_entity = *host_entity_pointer;

             LocalGridLocalFunction sub_loc_value = sub_func.localFunction( sub_entity );
             FineLocalFunction host_loc_value = host_func.localFunction( host_entity );

             const unsigned int numBaseFunctions = sub_loc_value.baseFunctionSet().numBaseFunctions();
             for( unsigned int i = 0; i < numBaseFunctions; ++i )
               {
                 host_loc_value[ i ] = sub_loc_value[ i ];
               }

          }

    }




  // dummy implementation of "operator()"
  // 'w' = effect of the discrete operator on 'u'
  template< class CoarseDiscreteFunctionImp, class MacroMicroGridSpecifierImp, class FineDiscreteFunctionImp, class DiffusionImp >
  void DiscreteEllipticMsFEMOperator< CoarseDiscreteFunctionImp, 
                                      MacroMicroGridSpecifierImp,
                                      FineDiscreteFunctionImp, 
                                      DiffusionImp > :: operator() ( const CoarseDiscreteFunction &u, CoarseDiscreteFunction &w ) const 
  {

    std :: cout << "the ()-operator of the DiscreteEllipticMsFEMOperator class is not yet implemented and still a dummy." << std :: endl;
    std :: abort();

  }


  template< class CoarseDiscreteFunctionImp, class MacroMicroGridSpecifierImp, class FineDiscreteFunctionImp, class DiffusionImp >
  template< class MatrixType >
  void DiscreteEllipticMsFEMOperator< CoarseDiscreteFunctionImp,
                                      MacroMicroGridSpecifierImp,
                                      FineDiscreteFunctionImp,
                                      DiffusionImp > :: assemble_matrix ( MatrixType &global_matrix ) const
  {
   
    // the local problem:
    // Let 'T' denote a coarse grid element and
    // let 'U(T)' denote the environment of 'T' that corresponds with the subgrid.
    
    // if Petrov-Galerkin-MsFEM
    #ifdef PGF
    std :: cout << "Assembling Petrov-Galerkin-MsFEM Matrix." << std :: endl;
    #else
    std :: cout << "Assembling MsFEM Matrix." << std :: endl;
    #endif



    //! der braucht einen macro-space (der wird auch als subgrid-space generiert) und den total globalen Feinskalen-Raum

    typedef typename MatrixType::LocalMatrixType LocalMatrix;

    global_matrix.reserve();
    global_matrix.clear();

    std::vector< typename CoarseBaseFunctionSet::JacobianRangeType > gradient_Phi( coarseDiscreteFunctionSpace_.mapper().maxNumDofs() );
   
    int local_problem_id = 0;

    // Coarse Entity Iterator 
    const CoarseIterator coarse_grid_end = coarseDiscreteFunctionSpace_.end();
    for( CoarseIterator coarse_grid_it = coarseDiscreteFunctionSpace_.begin(); coarse_grid_it != coarse_grid_end; ++coarse_grid_it )
    {

      // the coarse grid element T:
      const CoarseEntity &coarse_grid_entity = *coarse_grid_it;
      const CoarseGeometry &coarse_grid_geometry = coarse_grid_entity.geometry();
      assert( coarse_grid_entity.partitionType() == InteriorEntity );
    
      int global_index_entity = coarseDiscreteFunctionSpace_.gridPart().indexSet().index( coarse_grid_entity );
      
      LocalMatrix local_matrix = global_matrix.localMatrix( coarse_grid_entity, coarse_grid_entity );

      const CoarseBaseFunctionSet &coarse_grid_baseSet = local_matrix.domainBaseFunctionSet();
      const unsigned int numMacroBaseFunctions = coarse_grid_baseSet.numBaseFunctions();


#if 0
if ( global_index_entity == 0 )
{
std :: cout << "Im Assembler." << std :: endl;
std :: cout << "coarse_grid_it->geometry().corner(0) = " << coarse_grid_it->geometry().corner(0) << std :: endl;
std :: cout << "coarse_grid_it->geometry().corner(1) = " << coarse_grid_it->geometry().corner(1) << std :: endl;
std :: cout << "coarse_grid_it->geometry().corner(2) = " << coarse_grid_it->geometry().corner(2) << std :: endl;
}
#endif

#if 1

      // the sub grid U(T) that belongs to the coarse_grid_entity T
      SubGridType& sub_grid_U_T = subgrid_list_.getSubGrid( global_index_entity );
      SubGridPart subGridPart( sub_grid_U_T );
      
      LocalDiscreteFunctionSpace localDiscreteFunctionSpace( subGridPart );
      
      LocalDiscreteFunction local_problem_solution_e0( "Local problem Solution e_0", localDiscreteFunctionSpace );
      local_problem_solution_e0.clear();

      LocalDiscreteFunction local_problem_solution_e1( "Local problem Solution e_1", localDiscreteFunctionSpace );
      local_problem_solution_e1.clear();

      // --------- load local solutions -------

      char location_lps[50];
      sprintf( location_lps, "/_localProblemSolutions_%d", global_index_entity );
      std::string location_lps_s( location_lps );

      std :: string local_solution_location;

      // the file/place, where we saved the solutions of the cell problems
      local_solution_location = path_ + location_lps_s;

      bool reader_is_open = false;
      // reader for the cell problem data file:
      DiscreteFunctionReader discrete_function_reader( (local_solution_location).c_str() );
      reader_is_open = discrete_function_reader.open();

      if (reader_is_open)
        { discrete_function_reader.read( 0, local_problem_solution_e0 ); }

      if (reader_is_open)
        { discrete_function_reader.read( 1, local_problem_solution_e1 ); }

#endif

      // 1 point quadrature!! We only need the gradient of the base function,
      // which is constant on the whole entity.
      CoarseQuadrature one_point_quadrature( coarse_grid_entity, 0 );

      // the barycenter of the macro_grid_entity
      const typename CoarseQuadrature::CoordinateType &local_coarse_point 
           = one_point_quadrature.point( 0 /*=quadraturePoint*/ );
      DomainType coarse_entity_barycenter = coarse_grid_geometry.global( local_coarse_point );

      // transposed of the the inverse jacobian
      const FieldMatrix< double, dimension, dimension > &inverse_jac
          = coarse_grid_geometry.jacobianInverseTransposed( local_coarse_point );

      for( unsigned int i = 0; i < numMacroBaseFunctions; ++i )
        {
          // jacobian of the base functions, with respect to the reference element
          typename CoarseBaseFunctionSet::JacobianRangeType gradient_Phi_ref_element;
          coarse_grid_baseSet.jacobian( i, one_point_quadrature[ 0 ], gradient_Phi_ref_element );

          // multiply it with transpose of jacobian inverse to obtain the jacobian with respect to the real entity
          inverse_jac.mv( gradient_Phi_ref_element[ 0 ], gradient_Phi[ i ][ 0 ] );
        }

      for( unsigned int i = 0; i < numMacroBaseFunctions; ++i )
        {

          for( unsigned int j = 0; j < numMacroBaseFunctions; ++j )
           {

            RangeType local_integral = 0.0;
   
            // iterator for the micro grid ( grid for the reference element T_0 )
            const LocalGridIterator local_grid_end = localDiscreteFunctionSpace.end();
            for( LocalGridIterator local_grid_it = localDiscreteFunctionSpace.begin(); local_grid_it != local_grid_end; ++local_grid_it )
              {
		
                const LocalGridEntity &local_grid_entity = *local_grid_it;
		
		// check if "local_grid_entity" (which is an entity of U(T)) is in T:
		// -------------------------------------------------------------------

                FineEntityPointer fine_entity_pointer_1 = coarseDiscreteFunctionSpace_.grid().template getHostEntity<0>( coarse_grid_entity );
                FineEntityPointer fine_entity_pointer_2 = localDiscreteFunctionSpace.grid().template getHostEntity<0>( local_grid_entity );
		int coarse_level = coarseDiscreteFunctionSpace_.gridPart().grid().maxLevel();
		int fine_level = fineDiscreteFunctionSpace_.gridPart().grid().maxLevel();

                for (int lev = 0; lev < ( fine_level - coarse_level) ; ++lev)
                       fine_entity_pointer_2 = fine_entity_pointer_2->father();
		
		 bool entities_identical = true;
                int number_of_nodes = (*fine_entity_pointer_1).template count<2>();
                for ( int k = 0; k < number_of_nodes; k += 1 )
                  {
		     if ( !(fine_entity_pointer_1->geometry().corner(k) == fine_entity_pointer_2->geometry().corner(k)) )
                      { entities_identical = false; }
		  }

		 if ( entities_identical == false )
		  {
                   // std :: cout << "fine_entity_pointer_1->geometry().corner(0) = " << fine_entity_pointer_1->geometry().corner(0) << std :: endl;
                   // std :: cout << "fine_entity_pointer_1->geometry().corner(1) = " << fine_entity_pointer_1->geometry().corner(1) << std :: endl;
                   // std :: cout << "fine_entity_pointer_1->geometry().corner(2) = " << fine_entity_pointer_1->geometry().corner(2) << std :: endl;
                   // std :: cout << "fine_entity_pointer_2->geometry().corner(0) = " << fine_entity_pointer_2->geometry().corner(0) << std :: endl;
                   // std :: cout << "fine_entity_pointer_2->geometry().corner(1) = " << fine_entity_pointer_2->geometry().corner(1) << std :: endl;
                   // std :: cout << "fine_entity_pointer_2->geometry().corner(2) = " << fine_entity_pointer_2->geometry().corner(2) << std :: endl << std :: endl;
		   continue; 
		  }
       

		 // if ( !( fine_entity_pointer_1 == fine_entity_pointer_2 ) )
		 // { continue; }

		// -------------------------------------------------------------------
		
                const LocalGridGeometry &local_grid_geometry = local_grid_entity.geometry();
                assert( local_grid_entity.partitionType() == InteriorEntity );

                // higher order quadrature, since A^{\epsilon} is highly variable
                LocalGridQuadrature local_grid_quadrature( local_grid_entity, 2*localDiscreteFunctionSpace.order()+2 );
                const size_t numQuadraturePoints = local_grid_quadrature.nop();
		
                for( size_t localQuadraturePoint = 0; localQuadraturePoint < numQuadraturePoints; ++localQuadraturePoint )
                 {
                    // local (barycentric) coordinates (with respect to entity)
                    const typename LocalGridQuadrature::CoordinateType &local_subgrid_point = local_grid_quadrature.point( localQuadraturePoint );
		    
		    DomainType global_point_in_U_T = local_grid_geometry.global( local_subgrid_point );
		    
                    const double weight_local_quadrature 
                       = local_grid_quadrature.weight(  localQuadraturePoint ) * local_grid_geometry.integrationElement( local_subgrid_point );
		       
                    LocalGridLocalFunction localized_local_problem_solution_e0 = local_problem_solution_e0.localFunction( local_grid_entity );
                    LocalGridLocalFunction localized_local_problem_solution_e1 = local_problem_solution_e1.localFunction( local_grid_entity );

                    // grad coorector for e_0 and e_1
                    typename LocalGridBaseFunctionSet::JacobianRangeType grad_loc_sol_e0, grad_loc_sol_e1;
                    localized_local_problem_solution_e0.jacobian( local_grid_quadrature[ localQuadraturePoint ], grad_loc_sol_e0 );
                    localized_local_problem_solution_e1.jacobian( local_grid_quadrature[ localQuadraturePoint ], grad_loc_sol_e1 );
		    
		    // ∇ Phi_H + ∇ Q( Phi_H ) = ∇ Phi_H + ∂_x1 Phi_H Q( e_1 ) + ∂_x2 Phi_H Q( e_2 )
                    JacobianRangeType direction_of_diffusion( 0.0 );
                    for( int k = 0; k < dimension; ++k )
                      {
                        direction_of_diffusion[ 0 ][ k ] += gradient_Phi[ i ][ 0 ][ 0 ] * grad_loc_sol_e0[ 0 ][ k ];
                        direction_of_diffusion[ 0 ][ k ] += gradient_Phi[ i ][ 0 ][ 1 ] * grad_loc_sol_e1[ 0 ][ k ];
                        direction_of_diffusion[ 0 ][ k ] += gradient_Phi[ i ][ 0 ][ k ];
                      }
                      
                    JacobianRangeType diffusive_flux( 0.0 );
                    diffusion_operator_.diffusiveFlux( global_point_in_U_T, direction_of_diffusion, diffusive_flux );
                    
		    // if not Petrov-Galerkin:
                    #ifndef PGF
                    JacobianRangeType reconstruction_grad_phi_j( 0.0 );
                    for( int k = 0; k < dimension; ++k )
                      {
                        reconstruction_grad_phi_j[ 0 ][ k ] += gradient_Phi[ j ][ 0 ][ 0 ] * grad_loc_sol_e0[ 0 ][ k ];
                        reconstruction_grad_phi_j[ 0 ][ k ] += gradient_Phi[ j ][ 0 ][ 1 ] * grad_loc_sol_e1[ 0 ][ k ];
                        reconstruction_grad_phi_j[ 0 ][ k ] += gradient_Phi[ j ][ 0 ][ k ];
                      }

                    local_integral += weight_local_quadrature * ( diffusive_flux[ 0 ] *  reconstruction_grad_phi_j[ 0 ]);
                    #else
                    local_integral += weight_local_quadrature * ( diffusive_flux[ 0 ] * gradient_Phi[ j ][ 0 ]);
                    #endif

		  }
	       }

            // add entries
            local_matrix.add( j, i, local_integral );
           }

        }

    }
  
    
    // discrete_function_reader.close();

    // boundary treatment

    for( CoarseIterator coarse_grid_it = coarseDiscreteFunctionSpace_.begin(); coarse_grid_it != coarse_grid_end; ++coarse_grid_it )
    {
      
      const CoarseEntity &coarse_grid_entity = *coarse_grid_it;
      FineEntityPointer fine_entity_pointer = coarseDiscreteFunctionSpace_.grid().template getHostEntity<0>( coarse_grid_entity );
      
      const FineEntity& fine_entity = *fine_entity_pointer;
      
      LocalMatrix local_matrix = global_matrix.localMatrix( coarse_grid_entity, coarse_grid_entity );
      
      const CoarseLagrangePointSet &lagrangePointSet = coarseDiscreteFunctionSpace_.lagrangePointSet( coarse_grid_entity );
      
      const FineIntersectionIterator iend = fineDiscreteFunctionSpace_.gridPart().iend( fine_entity );
      for( FineIntersectionIterator iit = fineDiscreteFunctionSpace_.gridPart().ibegin( fine_entity ); iit != iend; ++iit )
        {
	  
           if ( iit->neighbor() ) //if there is a neighbor entity
            {
              // check if the neighbor entity is in the subgrid
              const FineEntityPointer neighborFineEntityPointer = iit->outside();
              const FineEntity& neighborFineEntity = *neighborFineEntityPointer;
              if ( coarseDiscreteFunctionSpace_.grid().template contains<0>( neighborFineEntity ) )
               {
                 continue;
               }

            }
            
           const int face = (*iit).indexInInside();
           const CoarseFaceDofIterator fdend = lagrangePointSet.template endSubEntity< 1 >( face );
           for( CoarseFaceDofIterator fdit = lagrangePointSet.template beginSubEntity< 1 >( face ); fdit != fdend; ++fdit )
              local_matrix.unitRow( *fdit );	
        }
  
    }

  }


//! ------------------------------------------------------------------------------------------------
//! ------------------------------------------------------------------------------------------------


}

#endif // #ifndef DiscreteElliptic_HH
