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
  template< class CoarseDiscreteFunctionImp, class FineDiscreteFunctionImp, class DiffusionImp >
  class DiscreteEllipticMsFEMOperator
  : public Operator< typename CoarseDiscreteFunctionImp::RangeFieldType,
                     typename CoarseDiscreteFunctionImp::RangeFieldType,
		      CoarseDiscreteFunctionImp, CoarseDiscreteFunctionImp >
  {
    typedef DiscreteEllipticMsFEMOperator< CoarseDiscreteFunctionImp, FineDiscreteFunctionImp, DiffusionImp > This;

  public:
    
    typedef CoarseDiscreteFunctionImp CoarseDiscreteFunction;
    typedef FineDiscreteFunctionImp FineDiscreteFunction;
   
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
    typedef SubGridList< FineDiscreteFunction, SubGridType > SubGridListType;
    
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

//!-----------------------------------------------------------------------------------------

  public:

    DiscreteEllipticMsFEMOperator( const CoarseDiscreteFunctionSpace &coarseDiscreteFunctionSpace,
                                   const FineDiscreteFunctionSpace &fineDiscreteFunctionSpace,
                                   // number of layers per coarse grid entity T:  U(T) is created by enrichting T with n(T)-layers:
                                   const std :: vector < int >& number_of_layers,
                                   const DiffusionModel &diffusion_op,
                                   std :: ofstream& data_file,
                                   std :: string path = ""  )
    : coarseDiscreteFunctionSpace_( coarseDiscreteFunctionSpace ),
      fineDiscreteFunctionSpace_( fineDiscreteFunctionSpace ),
      number_of_layers_( number_of_layers ),
      diffusion_operator_( diffusion_op ),
      data_file_( &data_file ),
      path_( path )
    {
      bool silence = false;

      const int coarse_level = coarseDiscreteFunctionSpace_.gridPart().grid().maxLevel();

      subgrid_list_ = new SubGridListType( fineDiscreteFunctionSpace_ , number_of_layers_, coarse_level , silence );

      SubGridListType subgrid_list ( *subgrid_list_ );

      //SubGridListType sl_( fineDiscreteFunctionSpace_ , number_of_layers_, coarse_level , silence );      
      //sl(sl_);
      //! Auslagern!!!!!!
      std :: string local_path =  path_ + "/local_problems/";
      MsFEMLocalProblemSolverType loc_prob_solver( fineDiscreteFunctionSpace_, subgrid_list, diffusion_operator_, data_file, local_path );
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
    
    const FineDiscreteFunctionSpace &fineDiscreteFunctionSpace_;
    const CoarseDiscreteFunctionSpace &coarseDiscreteFunctionSpace_;
    const std :: vector < int >& number_of_layers_;
    const DiffusionModel &diffusion_operator_;

    const SubGridListType* subgrid_list_;
    
    // data file for saving information
    std :: ofstream *data_file_;

    // path where to save the data output
    std :: string path_;

  };



  // dummy implementation of "operator()"
  // 'w' = effect of the discrete operator on 'u'
  template< class CoarseDiscreteFunctionImp, class FineDiscreteFunctionImp, class DiffusionImp >
  void DiscreteEllipticMsFEMOperator< CoarseDiscreteFunctionImp, 
                                      FineDiscreteFunctionImp, 
				       DiffusionImp > :: operator() ( const CoarseDiscreteFunction &u, CoarseDiscreteFunction &w ) const 
  {

    std :: cout << "the ()-operator of the DiscreteEllipticMsFEMOperator class is not yet implemented and still a dummy." << std :: endl;
    std :: abort();

  }


  template< class CoarseDiscreteFunctionImp, class FineDiscreteFunctionImp, class DiffusionImp >
  template< class MatrixType >
  void DiscreteEllipticMsFEMOperator< CoarseDiscreteFunctionImp, 
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

    std :: string local_solution_location;
    
    //! der braucht einen macro-space (der wird auch als subgrid-space generiert) und den total globalen Feinskalen-Raum
    
    // the file/place, where we saved the solutions of the cell problems
    local_solution_location = path_ + "/local_problems/_localProblemSolutions";

    bool reader_is_open = false;
    // reader for the cell problem data file:
    DiscreteFunctionReader discrete_function_reader( (local_solution_location).c_str() );
    reader_is_open = discrete_function_reader.open();

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
    
      int global_index_entity =  coarseDiscreteFunctionSpace_.gridPart().indexSet().index( coarse_grid_entity );
      
      LocalMatrix local_matrix = global_matrix.localMatrix( coarse_grid_entity, coarse_grid_entity );

      const CoarseBaseFunctionSet &coarse_grid_baseSet = local_matrix.domainBaseFunctionSet();
      const unsigned int numMacroBaseFunctions = coarse_grid_baseSet.numBaseFunctions();
      
#if 1
      SubGridListType subgrid_list ( *subgrid_list_ );
#if 0 
      // the sub grid U(T) that belongs to the coarse_grid_entity T
      SubGridType& sub_grid_U_T = subgrid_list.getSubGrid( global_index_entity );
      SubGridPart subGridPart( sub_grid_U_T );
      
      LocalDiscreteFunctionSpace localDiscreteFunctionSpace( subGridPart );
      
      LocalDiscreteFunction local_problem_solution_e0( "Local problem Solution e_0", localDiscreteFunctionSpace );
      local_problem_solution_e0.clear();
      
      LocalDiscreteFunction local_problem_solution_e1( "Local problem Solution e_1", localDiscreteFunctionSpace );
      local_problem_solution_e1.clear();
     
      if (reader_is_open)
        { discrete_function_reader.read( local_problem_id, local_problem_solution_e0 ); }
      local_problem_id += 1;
      if (reader_is_open)
        { discrete_function_reader.read( local_problem_id, local_problem_solution_e1 ); }
#endif
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
	  
      // DiscreteFunction* corrector_Phi[discreteFunctionSpace_.mapper().maxNumDofs()];
	  
      //! TODO: get local grids, load local solutions     
	  
      for( unsigned int i = 0; i < numMacroBaseFunctions; ++i )
        {

          // jacobian of the base functions, with respect to the reference element
          typename CoarseBaseFunctionSet::JacobianRangeType gradient_Phi_ref_element;
          coarse_grid_baseSet.jacobian( i, one_point_quadrature[ 0 ], gradient_Phi_ref_element );

          // multiply it with transpose of jacobian inverse to obtain the jacobian with respect to the real entity
          inverse_jac.mv( gradient_Phi_ref_element[ 0 ], gradient_Phi[ i ][ 0 ] );

          // --------- load local solutions -------
          // corrector_Phi[i] = new DiscreteFunction( "Corrector Function of Phi", localDiscreteFunctionSpace_ );
          // corrector_Phi[i]->clear();
          //if (reader_is_open)
          //  { discrete_function_reader.read( cell_problem_id[ i ], *(corrector_Phi[ i ]) ); }


        }

      for( unsigned int i = 0; i < numMacroBaseFunctions; ++i )
        {

          for( unsigned int j = 0; j < numMacroBaseFunctions; ++j )
           {

            RangeType local_integral = 0.0;
    
#if 0

   
#if 0
      




#if 1
DomainType A_0_MSFEM(0.0);
#endif

            // iterator for the micro grid ( grid for the reference element T_0 )
            const Iterator micro_grid_end = localDiscreteFunctionSpace_.end();
            for( Iterator micro_grid_it = localDiscreteFunctionSpace_.begin(); micro_grid_it != micro_grid_end; ++micro_grid_it )
              {

                // remember:
                // |det(A)| \int_{T_0} (A^eps ○ F)(x) ( ∇\Phi_i(x_T) + (A^{-1})^T ∇( Q^eps(\Phi_i) ○ F )(x)) · ( ∇\Phi_j(x_T) + (A^{-1})^T ∇( Q^eps(\Phi_j) ○ F )(x))

#if 1
gradient_Phi[ i ][ 0 ][ 0 ] = 1.0;
gradient_Phi[ i ][ 0 ][ 1 ] = 0.0;
#endif

                const Entity &micro_grid_entity = *micro_grid_it;
                const Geometry &micro_grid_geometry = micro_grid_entity.geometry();
                assert( micro_grid_entity.partitionType() == InteriorEntity );

                // ( Q^eps(\Phi_i) ○ F ):
                typename DiscreteFunction::LocalFunctionType localized_corrector_i = corrector_Phi[ i ]->localFunction( micro_grid_entity );
                // ( Q^eps(\Phi_j) ○ F ):
                typename DiscreteFunction::LocalFunctionType localized_corrector_j = corrector_Phi[ j ]->localFunction( micro_grid_entity );

                // higher order quadrature, since A^{\epsilon} is highly variable
                Quadrature micro_grid_quadrature( micro_grid_entity, 2*localDiscreteFunctionSpace_.order()+2 );
                const size_t numQuadraturePoints = micro_grid_quadrature.nop();

                for( size_t microQuadraturePoint = 0; microQuadraturePoint < numQuadraturePoints; ++microQuadraturePoint )
                  {

                    // local (barycentric) coordinates (with respect to entity)
                    const typename Quadrature::CoordinateType &local_micro_point = micro_grid_quadrature.point( microQuadraturePoint );

                    DomainType global_point_in_T_0 = micro_grid_geometry.global( local_micro_point );

                    const double weight_micro_quadrature = micro_grid_quadrature.weight(  microQuadraturePoint ) * micro_grid_geometry.integrationElement( local_micro_point );

                    // ∇( Q^eps(\Phi_i) ○ F )   and   ∇( Q^eps(\Phi_j) ○ F )
                    JacobianRangeType grad_corrector_i, grad_corrector_j;
                    localized_corrector_i.jacobian( micro_grid_quadrature[ microQuadraturePoint ], grad_corrector_i );
                    localized_corrector_j.jacobian( micro_grid_quadrature[ microQuadraturePoint ], grad_corrector_j );

#if 1
                    // global point in the reference element T_0
                    DomainType global_point = micro_grid_geometry.global( local_micro_point );

                    // 'F(x)', i.e. F ( global point in the reference element T_0 )
                    // (the transformation of the global point in T_0 to its position in T)
                    DomainType global_point_transformed(0.0);

                    for( int k = 0; k < dimension; ++k )
                      for( int l = 0; l < dimension; ++l )
                        global_point_transformed[ k ] += ( val_A[ k ][ l ] * global_point_in_T_0[ l ] );

                    global_point_transformed += corner_0_of_T;

                    // F(x) = Ax + a_0, F : T_0 -> T is given by
                    // A_11 = a_1(1) - a_0(1)     A_12 = a_2(1) - a_0(1)
                    // A_21 = a_1(2) - a_0(2)     A_22 = a_2(2) - a_0(2)

#endif
                    // direction of the diffusion: ( ∇\Phi_i(x_T) + (A^{-1})^T ∇( Q^eps(\Phi_i) ○ F )(x))
                    JacobianRangeType direction_of_diffusion( 0.0 );
                    for( int k = 0; k < dimension; ++k )
                      {
                        for( int l = 0; l < dimension; ++l )
                          {
                           direction_of_diffusion[ 0 ][ k ] += val_A_inverse_transposed[ k ][ l ] * grad_corrector_i[ 0 ][ l ];
                          }

                        direction_of_diffusion[ 0 ][ k ] += gradient_Phi[ i ][ 0 ][ k ];

                      }

                    JacobianRangeType diffusion_in_gradient_Phi_reconstructed;
                    diffusion_operator_.diffusiveFlux( global_point_transformed,
                                                       direction_of_diffusion, diffusion_in_gradient_Phi_reconstructed );


                    // if test function reconstruction
                    #ifndef PGF


                    // ( ∇\Phi_j(x_T) + (A^{-1})^T ∇( Q^eps(\Phi_j) ○ F )(x)):
                    JacobianRangeType grad_reconstruction_Phi_j( 0.0 );
                    for( int k = 0; k < dimension; ++k )
                      {
                        for( int l = 0; l < dimension; ++l )
                          { grad_reconstruction_Phi_j[ 0 ][ k ] += val_A_inverse_transposed[ k ][ l ] * grad_corrector_j[ 0 ][ l ]; }
                        grad_reconstruction_Phi_j[ 0 ][ k ] += gradient_Phi[ j ][ 0 ][ k ];
                      }

                    local_integral += weight_micro_quadrature * ( diffusion_in_gradient_Phi_reconstructed[ 0 ] * grad_reconstruction_Phi_j[ 0 ]);
                    #else
                    local_integral += weight_micro_quadrature * ( diffusion_in_gradient_Phi_reconstructed[ 0 ] * gradient_Phi[ j ][ 0 ]);

#if 1
A_0_MSFEM[ 0 ] +=  2.0 * weight_micro_quadrature * diffusion_in_gradient_Phi_reconstructed[ 0 ][ 0 ];
A_0_MSFEM[ 1 ] +=  2.0 * weight_micro_quadrature * diffusion_in_gradient_Phi_reconstructed[ 0 ][ 1 ];
#endif

                    #endif

                  }
              }

#if 1
std :: cout << "A_0_MSFEM[ 0 ] = " << A_0_MSFEM[ 0 ] << std :: endl;
std :: cout << "A_0_MSFEM[ 1 ] = " << A_0_MSFEM[ 1 ] << std :: endl;
#endif

#endif
#endif
            // add |det(A)|*\int_{T_0} ...
            local_matrix.add( j, i, /*! .... values ... */ local_integral );
           }

        }

      // delete?
      //delete[] corrector_Phi;

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
