#ifndef DiscreteEllipticMSFEMOperator_HH
#define DiscreteEllipticMSFEMOperator_HH

#include <dune/common/fmatrix.hh>

#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/multiscale/tools/solver/MsFEM/msfem_localproblems/localproblemsolver.hh>

#include <dune/fem/operator/2order/lagrangematrixsetup.hh>

namespace Dune
{

  // Imp stands for Implementation
  template< class DiscreteFunctionImp, class DiffusionImp >
  class DiscreteEllipticMsFEMOperator
  : public Operator< typename DiscreteFunctionImp::RangeFieldType, typename DiscreteFunctionImp::RangeFieldType, DiscreteFunctionImp, DiscreteFunctionImp >
  {
    typedef DiscreteEllipticMsFEMOperator< DiscreteFunctionImp, DiffusionImp > This;

  public:

    typedef DiscreteFunctionImp DiscreteFunction;
   
    typedef DiffusionImp DiffusionModel;

    typedef typename DiscreteFunction::DiscreteFunctionSpaceType DiscreteFunctionSpace;

    typedef typename DiscreteFunctionSpace::GridPartType GridPart;
    typedef typename DiscreteFunctionSpace::GridType GridType;

    typedef typename DiscreteFunctionSpace::RangeFieldType RangeFieldType;
    typedef typename DiscreteFunctionSpace::DomainType DomainType;
    typedef typename DiscreteFunctionSpace::RangeType RangeType;
    typedef typename DiscreteFunctionSpace::JacobianRangeType
      JacobianRangeType;

    typedef MsFEMLocalProblemSolver< DiscreteFunction, DiffusionModel > MsFEMLocalProblemSolverType;

  protected:
    static const int dimension = GridPart::GridType::dimension;
    static const int polynomialOrder = DiscreteFunctionSpace::polynomialOrder;

    typedef typename DiscreteFunction::LocalFunctionType LocalFunction;

    typedef typename DiscreteFunctionSpace::BaseFunctionSetType BaseFunctionSet;
    typedef typename DiscreteFunctionSpace::LagrangePointSetType LagrangePointSet;
    typedef typename LagrangePointSet::template Codim< 1 >::SubEntityIteratorType FaceDofIterator;

    typedef typename DiscreteFunctionSpace::IteratorType Iterator;
    typedef typename Iterator::Entity Entity;
    typedef typename Entity::Geometry Geometry;

    typedef typename GridPart::IntersectionIteratorType IntersectionIterator;
    typedef typename IntersectionIterator::Intersection Intersection;

    typedef typename GridType :: template Codim< 0 > :: template  Partition< All_Partition > :: LevelIterator LevelEntityIterator;

    typedef CachingQuadrature< GridPart, 0 > Quadrature;

  public:
    
    DiscreteEllipticMsFEMOperator( const DiscreteFunctionSpace &discreteFunctionSpace,
				    const int coarse_grid_level, // Refinment Level for the coarse grid 
                                   DiffusionModel &diffusion_op,
                                   std :: string &filename = "default") //!NOTE: default ueberarbeiten!
    : discreteFunctionSpace_( discreteFunctionSpace ),
      coarse_grid_level_( coarse_grid_level ),
      diffusion_operator_( diffusion_op )
     //filename(filename)
    { /* typedef MsFEMLocalProblemSolver< DiscreteFunction, DiffusionModel > MsFEMLocalProblemSolverType;*/ }


  private:
    DiscreteEllipticMsFEMOperator ( const This & );

  public:

    // dummy operator
    virtual void
    operator() ( const DiscreteFunction &u, DiscreteFunction &w ) const;


    template< class MatrixType >
    void assemble_matrix ( MatrixType &global_matrix ) const;

  private:
    const DiscreteFunctionSpace &discreteFunctionSpace_;
    DiffusionModel &diffusion_operator_;

    const int coarse_grid_level_;
    
    // name of data file, e.g. required if we want to use the saved solutions of the cell problems
    std :: string *filename_;

  };



  // dummy implementation of "operator()"
  // 'w' = effect of the discrete operator on 'u'
  template< class DiscreteFunctionImp, class DiffusionImp >
  void DiscreteEllipticMsFEMOperator< DiscreteFunctionImp, DiffusionImp >::operator() ( const DiscreteFunction &u, DiscreteFunction &w ) const 
  {

    std :: cout << "the ()-operator of the DiscreteEllipticMsFEMOperator class is not yet implemented and still a dummy." << std :: endl;
    std :: abort();

  }


  template< class DiscreteFunctionImp, class DiffusionImp >
  template< class MatrixType >
  void DiscreteEllipticMsFEMOperator< DiscreteFunctionImp, DiffusionImp >::assemble_matrix ( MatrixType &global_matrix ) const
  {
#if 0
    
    //! der braucht einen macro-space (der wird auch als subgrid-space generiert) und den total globalen Feinskalen-Raum
    
   
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

    // the file/place, where we saved the solutions of the cell problems
    local_solution_location = "data/MsFEM/"+(*filename_)+"/local_problems/_localProblemSolutions_baseSet";


    bool reader_is_open = false;
    // reader for the cell problem data file:
    DiscreteFunctionReader discrete_function_reader( (local_solution_location).c_str() );
    reader_is_open = discrete_function_reader.open();

    typedef typename MatrixType::LocalMatrixType LocalMatrix;

    global_matrix.reserve();
    global_matrix.clear();

    std::vector< typename BaseFunctionSet::JacobianRangeType > gradient_Phi( discreteFunctionSpace_.mapper().maxNumDofs() );
    
//! assemble local grids, load local solutions
    
    //! LevelEntityIterator
    
    const Iterator macro_grid_end = discreteFunctionSpace_.end();
    for( Iterator macro_grid_it = discreteFunctionSpace_.begin(); macro_grid_it != macro_grid_end; ++macro_grid_it )
    {

      // the coarse grid element T:
      const Entity &macro_grid_entity = *macro_grid_it;
      const Geometry &macro_grid_geometry = macro_grid_entity.geometry();
      assert( macro_grid_entity.partitionType() == InteriorEntity );

      LocalMatrix local_matrix = global_matrix.localMatrix( macro_grid_entity, macro_grid_entity );

      const BaseFunctionSet &macro_grid_baseSet = local_matrix.domainBaseFunctionSet();
      const unsigned int numMacroBaseFunctions = macro_grid_baseSet.numBaseFunctions();

#if 0
      
      // 1 point quadrature!! That is how we compute and save the cell problems.
      // If you want to use a higher order quadrature, you also need to change the computation of the cell problems!
      Quadrature one_point_quadrature( macro_grid_entity, 0 );

      // the barycenter of the macro_grid_entity
      const typename Quadrature::CoordinateType &local_macro_point = one_point_quadrature.point( 0 /*=quadraturePoint*/ );
      DomainType macro_entity_barycenter = macro_grid_geometry.global( local_macro_point );

      const double macro_entity_volume = one_point_quadrature.weight( 0 /*=quadraturePoint*/ ) * macro_grid_geometry.integrationElement( local_macro_point );

      // transposed of the the inverse jacobian
      const FieldMatrix< double, dimension, dimension > &inverse_jac
          = macro_grid_geometry.jacobianInverseTransposed( local_macro_point );

      int cell_problem_id [ numMacroBaseFunctions ];

      DiscreteFunction* corrector_Phi[discreteFunctionSpace_.mapper().maxNumDofs()];

      for( unsigned int i = 0; i < numMacroBaseFunctions; ++i )
        {

          // get number of cell problem from entity and number of base function
          cell_problem_id[i] = lp_num_manager_.get_number_of_local_problem( macro_grid_it, i );

          // jacobian of the base functions, with respect to the reference element
          typename BaseFunctionSet::JacobianRangeType gradient_Phi_ref_element;
          macro_grid_baseSet.jacobian( i, one_point_quadrature[ 0 ], gradient_Phi_ref_element );

          // multiply it with transpose of jacobian inverse to obtain the jacobian with respect to the real entity
          inverse_jac.mv( gradient_Phi_ref_element[ 0 ], gradient_Phi[ i ][ 0 ] );

          // ( Q^eps(\Phi_i) ○ F ):
          corrector_Phi[i] = new DiscreteFunction( "Corrector Function of Phi", localDiscreteFunctionSpace_ );
          corrector_Phi[i]->clear();
          #ifdef AD_HOC_COMPUTATION
          MsFEMLocalProblemSolverType cell_problem_solver( localDiscreteFunctionSpace_, diffusion_operator_ );
          cell_problem_solver.template solvelocalproblem<JacobianRangeType>
                ( gradient_Phi[ i ], macro_entity_barycenter, *(corrector_Phi[ i ]) );
          #else
          if (reader_is_open)
            { discrete_function_reader.read( cell_problem_id[ i ], *(corrector_Phi[ i ]) ); }
          #endif

        }

      for( unsigned int i = 0; i < numMacroBaseFunctions; ++i )
        {

          for( unsigned int j = 0; j < numMacroBaseFunctions; ++j )
           {

            RangeType local_integral = 0.0;


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

            // add |det(A)|*\int_{T_0} ...
            local_matrix.add( j, i, abs_det_A * local_integral );
           }

        }

      // delete?
      //delete[] corrector_Phi;
#endif
    }

    //discrete_function_reader.close();

    // boundary treatment
    const GridPart &gridPart = discreteFunctionSpace_.gridPart();
    for( Iterator it = discreteFunctionSpace_.begin(); it != macro_grid_end; ++it )
    {
      const Entity &entity = *it;
      if( !entity.hasBoundaryIntersections() )
        continue;

      LocalMatrix local_matrix = global_matrix.localMatrix( entity, entity );

      const LagrangePointSet &lagrangePointSet = discreteFunctionSpace_.lagrangePointSet( entity );

      const IntersectionIterator iend = gridPart.iend( entity );
      for( IntersectionIterator iit = gridPart.ibegin( entity ); iit != iend; ++iit )
      {
        const Intersection &intersection = *iit;
        if( !intersection.boundary() )
          continue;

        const int face = intersection.indexInInside();
        const FaceDofIterator fdend = lagrangePointSet.template endSubEntity< 1 >( face );
        for( FaceDofIterator fdit = lagrangePointSet.template beginSubEntity< 1 >( face ); fdit != fdend; ++fdit )
          local_matrix.unitRow( *fdit );
      }
    }
#endif
  }


//! ------------------------------------------------------------------------------------------------
//! ------------------------------------------------------------------------------------------------


}

#endif // #ifndef DiscreteElliptic_HH
