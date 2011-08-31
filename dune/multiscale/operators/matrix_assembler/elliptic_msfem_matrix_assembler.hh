#ifndef DiscreteEllipticMSFEMOperator_HH
#define DiscreteEllipticMSFEMOperator_HH

#include <dune/common/fmatrix.hh>

#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/multiscale/operators/msfem_localproblems/localproblemsolver.hh>

#include <dune/fem/operator/2order/lagrangematrixsetup.hh>

namespace Dune
{

  // Imp stands for Implementation
  template< class DiscreteFunctionImp, class DiffusionImp, class LocalProblemNumberingManagerImp >
  class DiscreteEllipticMsFEMOperator
  : public Operator< typename DiscreteFunctionImp::RangeFieldType, typename DiscreteFunctionImp::RangeFieldType, DiscreteFunctionImp, DiscreteFunctionImp >
  {
    typedef DiscreteEllipticMsFEMOperator< DiscreteFunctionImp, DiffusionImp, LocalProblemNumberingManagerImp > This;

  public:

    typedef DiscreteFunctionImp DiscreteFunction;
    typedef DiffusionImp DiffusionModel;
    typedef LocalProblemNumberingManagerImp LocalProblemNumberingManager;

    typedef typename DiscreteFunction::DiscreteFunctionSpaceType DiscreteFunctionSpace;

    typedef typename DiscreteFunctionSpace::GridPartType GridPart;
    typedef typename DiscreteFunctionSpace::GridType GridType;

    typedef typename DiscreteFunctionSpace::RangeFieldType RangeFieldType;
    typedef typename DiscreteFunctionSpace::DomainType DomainType;
    typedef typename DiscreteFunctionSpace::RangeType RangeType;
    typedef typename DiscreteFunctionSpace::JacobianRangeType
      JacobianRangeType;

    #ifdef AD_HOC_COMPUTATION
    typedef MsFEMLocalProblemSolver< DiscreteFunction, DiffusionModel > MsFEMLocalProblemSolverType;
    #endif

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

    typedef CachingQuadrature< GridPart, 0 > Quadrature;

  public:
    DiscreteEllipticMsFEMOperator( const DiscreteFunctionSpace &discreteFunctionSpace,
                                   const DiscreteFunctionSpace &localDiscreteFunctionSpace,
                                   DiffusionModel &diffusion_op,
                                   LocalProblemNumberingManager &cp_num_manager )
    : discreteFunctionSpace_( discreteFunctionSpace ),
      localDiscreteFunctionSpace_( localDiscreteFunctionSpace ),
      diffusion_operator_( diffusion_op ),
      filename_( NULL )
    {}

    DiscreteEllipticMsFEMOperator( const DiscreteFunctionSpace &discreteFunctionSpace,
                                   const DiscreteFunctionSpace &localDiscreteFunctionSpace,
                                   DiffusionModel &diffusion_op,
                                   LocalProblemNumberingManager &cp_num_manager,
                                   std :: string &filename )
    : discreteFunctionSpace_( discreteFunctionSpace ),
      localDiscreteFunctionSpace_( localDiscreteFunctionSpace ),
      diffusion_operator_( diffusion_op ),
      cp_num_manager_( cp_num_manager ),
      filename_( &filename )
    {}


  private:
    DiscreteEllipticMsFEMOperator ( const This & );

  public:

    // dummy operator
    virtual void
    operator() ( const DiscreteFunction &u, DiscreteFunction &w ) const;


    template< class MatrixType >
    void assemble_matrix ( MatrixType &global_matrix ) const;

    template< class MatrixType >
    void assemble_jacobian_matrix ( DiscreteFunction &old_macro_function, MatrixType &global_matrix ) const;


  private:
    const DiscreteFunctionSpace &discreteFunctionSpace_;
    const DiscreteFunctionSpace &localDiscreteFunctionSpace_;
    DiffusionModel &diffusion_operator_;

    // name of data file, e.g. required if we want to use the saved solutions of the cell problems
    std :: string *filename_;
    LocalProblemNumberingManager &cp_num_manager_;
  };



  // dummy implementation of "operator()"
  // 'w' = effect of the discrete operator on 'u'
  template< class DiscreteFunctionImp, class DiffusionImp, class LocalProblemNumberingManagerImp >
  void DiscreteEllipticMsFEMOperator< DiscreteFunctionImp, DiffusionImp, LocalProblemNumberingManagerImp >::operator() ( const DiscreteFunction &u, DiscreteFunction &w ) const 
  {

    std :: cout << "the ()-operator of the DiscreteEllipticMsFEMOperator class is not yet implemented and still a dummy." << std :: endl;
    std :: abort();

  }


  template< class DiscreteFunctionImp, class DiffusionImp, class LocalProblemNumberingManagerImp >
  template< class MatrixType >
  void DiscreteEllipticMsFEMOperator< DiscreteFunctionImp, DiffusionImp, LocalProblemNumberingManagerImp >::assemble_matrix ( MatrixType &global_matrix ) const
  {

    // if test function reconstruction
    #ifdef PGF
    std :: cout << "Assembling Petrov-Galerkin-MsFEM Matrix." << std :: endl;
    #else
    std :: cout << "Assembling MsFEM Matrix." << std :: endl;
    #endif

    std :: string cell_solution_location;

    // if we know the file, where we saved the solutions of the cell problems, use it
    if ( filename_ )
     {
      // place, where we saved the solutions of the cell problems
      cell_solution_location = "data/MsFEM/"+(*filename_)+"/local_problems/_localProblemSolutions_baseSet";
     }
    else
     {
      #ifndef AD_HOC_COMPUTATION
      std :: cout << "ERROR! No 'filename_' in class 'DiscreteEllipticMsFEMOperator', but no AD_HOC_COMPUTATION initialized. Therefore the location of the saved cell problems is not available. Please define AD_HOC_COMPUTATION (ad hoc computation of the cell problems) or pass a corresponding 'filename_'-variable!" << std :: endl;
      std :: abort();
      #endif
     }

    Problem::ModelProblemData model_info;
    const double delta = model_info.getDelta();
    const double epsilon_estimated = model_info.getEpsilonEstimated();


    bool reader_is_open = false;
    // reader for the cell problem data file:
    DiscreteFunctionReader discrete_function_reader( (cell_solution_location).c_str() );
    reader_is_open = discrete_function_reader.open();

    typedef typename MatrixType::LocalMatrixType LocalMatrix;

    global_matrix.reserve();
    global_matrix.clear();

    std::vector< typename BaseFunctionSet::JacobianRangeType > gradient_Phi( discreteFunctionSpace_.mapper().maxNumDofs() );

    const Iterator macro_grid_end = discreteFunctionSpace_.end();
    for( Iterator macro_grid_it = discreteFunctionSpace_.begin(); macro_grid_it != macro_grid_end; ++macro_grid_it )
    {

      const Entity &macro_grid_entity = *macro_grid_it;
      const Geometry &macro_grid_geometry = macro_grid_entity.geometry();
      assert( macro_grid_entity.partitionType() == InteriorEntity );

      LocalMatrix local_matrix = global_matrix.localMatrix( macro_grid_entity, macro_grid_entity );

      const BaseFunctionSet &macro_grid_baseSet = local_matrix.domainBaseFunctionSet();
      const unsigned int numMacroBaseFunctions = macro_grid_baseSet.numBaseFunctions();

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
          cell_problem_id[i] = cp_num_manager_.get_number_of_local_problem( macro_grid_it, i );

          // jacobian of the base functions, with respect to the reference element
          typename BaseFunctionSet::JacobianRangeType gradient_Phi_ref_element;
          macro_grid_baseSet.jacobian( i, one_point_quadrature[ 0 ], gradient_Phi_ref_element );

          // multiply it with transpose of jacobian inverse to obtain the jacobian with respect to the real entity
          inverse_jac.mv( gradient_Phi_ref_element[ 0 ], gradient_Phi[ i ][ 0 ] );

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

            RangeType fine_scale_average = 0.0;

            // nur checken ob der momentane Quadraturpunkt in der Zelle delta/epsilon_estimated*Y ist (0 bzw. 1 => Abschneidefunktion!)

            const Iterator micro_grid_end = localDiscreteFunctionSpace_.end();
            for( Iterator micro_grid_it = localDiscreteFunctionSpace_.begin(); micro_grid_it != micro_grid_end; ++micro_grid_it )
              {

                const Entity &micro_grid_entity = *micro_grid_it;
                const Geometry &micro_grid_geometry = micro_grid_entity.geometry();
                assert( micro_grid_entity.partitionType() == InteriorEntity );

                typename DiscreteFunction::LocalFunctionType localized_corrector_i = corrector_Phi[ i ]->localFunction( micro_grid_entity );
                typename DiscreteFunction::LocalFunctionType localized_corrector_j = corrector_Phi[ j ]->localFunction( micro_grid_entity );

                // higher order quadrature, since A^{\epsilon} is highly variable
                Quadrature micro_grid_quadrature( micro_grid_entity, 2*localDiscreteFunctionSpace_.order()+2 );
                const size_t numQuadraturePoints = micro_grid_quadrature.nop();

                for( size_t microQuadraturePoint = 0; microQuadraturePoint < numQuadraturePoints; ++microQuadraturePoint )
                  {

                    // local (barycentric) coordinates (with respect to entity)
                    const typename Quadrature::CoordinateType &local_micro_point = micro_grid_quadrature.point( microQuadraturePoint );

                    DomainType global_point_in_Y = micro_grid_geometry.global( local_micro_point );

                    const double weight_micro_quadrature = micro_grid_quadrature.weight(  microQuadraturePoint ) * micro_grid_geometry.integrationElement( local_micro_point );

                    JacobianRangeType grad_corrector_i, grad_corrector_j;
                    localized_corrector_i.jacobian( micro_grid_quadrature[ microQuadraturePoint ], grad_corrector_i );
                    localized_corrector_j.jacobian( micro_grid_quadrature[ microQuadraturePoint ], grad_corrector_j );

                    // x_T + (delta * y)
                    DomainType current_point_in_macro_grid;
                    for( int k = 0; k < dimension; ++k )
                      current_point_in_macro_grid[ k ] = macro_entity_barycenter[ k ] + (delta * global_point_in_Y[ k ]);

                    JacobianRangeType direction_of_diffusion;
                    for( int k = 0; k < dimension; ++k )
                       direction_of_diffusion[ 0 ][ k ] = gradient_Phi[ i ][ 0 ][ k ] + grad_corrector_i[ 0 ][ k ];


                    JacobianRangeType diffusion_in_gradient_Phi_reconstructed;
                    diffusion_operator_.diffusiveFlux( current_point_in_macro_grid,
                                                       direction_of_diffusion, diffusion_in_gradient_Phi_reconstructed );


                    double cutting_function = 1.0;
                    for( int k = 0; k < dimension; ++k )
                      {
                        // is the current quadrature point in the relevant cell?
                        if ( fabs(global_point_in_Y[ k ]) > (0.5*(epsilon_estimated/delta)) )
                           { cutting_function *= 0.0; }
                      }

                    // if test function reconstruction
                    #ifndef PGF
                    JacobianRangeType grad_reconstruction_Phi_j;
                    for( int k = 0; k < dimension; ++k )
                       grad_reconstruction_Phi_j[ 0 ][ k ] = gradient_Phi[ j ][ 0 ][ k ] + grad_corrector_j[ 0 ][ k ];

                    fine_scale_average += cutting_function * weight_micro_quadrature * ( diffusion_in_gradient_Phi_reconstructed[ 0 ] * grad_reconstruction_Phi_j[ 0 ]);
                    #else
                    fine_scale_average += cutting_function * weight_micro_quadrature * ( diffusion_in_gradient_Phi_reconstructed[ 0 ] * gradient_Phi[ j ][ 0 ]);
                    #endif

                  }
              }

            // add |T| * (delta/epsilon)^N \int_Y ...
            local_matrix.add( j, i,
                 pow( delta / epsilon_estimated, dimension ) * macro_entity_volume * fine_scale_average );
           }

        }

      // delete?
      //delete[] corrector_Phi;

    }

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
  }

//! ------------------------------------------------------------------------------------------------
//! ------------------------------------------------------------------------------------------------



}

#endif // #ifndef DiscreteElliptic_HH
