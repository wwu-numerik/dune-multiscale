#ifndef DiscreteElliptic_HH
#define DiscreteElliptic_HH

#include <dune/common/fmatrix.hh>

#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/operator/common/operator.hh>

namespace Dune
{

  // Imp stands for Implementation
  template< class DiscreteFunctionImp, class DiffusionImp, class ReactionImp >
  class DiscreteEllipticOperator
  : public Operator< typename DiscreteFunctionImp::RangeFieldType, typename DiscreteFunctionImp::RangeFieldType, DiscreteFunctionImp, DiscreteFunctionImp >
  {
    typedef DiscreteEllipticOperator< DiscreteFunctionImp, DiffusionImp, ReactionImp > This;

  public:
    typedef DiscreteFunctionImp DiscreteFunction;
    typedef DiffusionImp DiffusionModel;
    typedef ReactionImp Reaction;

    typedef typename DiscreteFunction::DiscreteFunctionSpaceType DiscreteFunctionSpace;

    typedef typename DiscreteFunctionSpace::GridPartType GridPart;
    typedef typename DiscreteFunctionSpace::RangeFieldType RangeFieldType;

    typedef typename DiscreteFunctionSpace::DomainType DomainType;
    typedef typename DiscreteFunctionSpace::RangeType RangeType;

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
    DiscreteEllipticOperator( const DiscreteFunctionSpace &discreteFunctionSpace, const DiffusionModel &diffusion_op )
    : discreteFunctionSpace_( discreteFunctionSpace ),
      diffusion_operator_( diffusion_op )
    {}

    DiscreteEllipticOperator( const DiscreteFunctionSpace &discreteFunctionSpace, const DiffusionModel &diffusion_op, Reaction &reaction_coefficient )
    : discreteFunctionSpace_( discreteFunctionSpace ),
      diffusion_operator_( diffusion_op ),
      reaction_coefficient_( &reaction_coefficient )
    {}

        
  private:
    DiscreteEllipticOperator ( const This & );

  public:

    // dummy operator
    virtual void
    operator() ( const DiscreteFunction &u, DiscreteFunction &w ) const;


    template< class MatrixType >
    void assemble_matrix ( MatrixType &global_matrix, bool boundary_treatment ) const;

    template< class MatrixType >
    void assemble_jacobian_matrix ( DiscreteFunction &disc_func, MatrixType &global_matrix, bool boundary_treatment ) const;


  private:
    const DiscreteFunctionSpace &discreteFunctionSpace_;
    const DiffusionModel &diffusion_operator_;
    Reaction *reaction_coefficient_;
  };


  // dummy implementation of "operator()"
  // 'w' = effect of the discrete operator on 'u'
  template< class DiscreteFunctionImp, class DiffusionImp, class ReactionImp >
  void DiscreteEllipticOperator< DiscreteFunctionImp, DiffusionImp, ReactionImp >::operator() ( const DiscreteFunction &u, DiscreteFunction &w ) const 
  {

    std :: cout << "the ()-operator of the DiscreteEllipticOperator class is not yet implemented and still a dummy." << std :: endl;
    std :: abort();

    // for the discrete Laplacian operator, we might use the following version:
#if 0

    w.clear();

    const Iterator end = discreteFunctionSpace_.end();
    for( Iterator it = discreteFunctionSpace_.begin(); it != end; ++it )
    {
      const Entity &entity = *it;
      const Geometry &geometry = entity.geometry();
      assert( entity.partitionType() == InteriorEntity );

      const LocalFunction uLocal = u.localFunction( entity );
      LocalFunction wLocal = w.localFunction( entity );
      
      const BaseFunctionSet &baseSet = wLocal.baseFunctionSet();
      const unsigned int numBaseFunctions = baseSet.numBaseFunctions();
            
      // for constant diffusion "2*discreteFunctionSpace_.order()" is sufficient, for the general case, it is better to use a higher order quadrature:
      Quadrature quadrature( entity, 2*discreteFunctionSpace_.order()+2 );
      const size_t numQuadraturePoints = quadrature.nop();
      for( size_t quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint )
      {
        // local (barycentric) coordinates (with respect to entity)
        const typename Quadrature::CoordinateType &local_point = quadrature.point( quadraturePoint );

        DomainType global_point = geometry.global( local_point );

        const double weight = quadrature.weight( quadraturePoint ) * geometry.integrationElement( local_point );
        const FieldMatrix< double, dimension, dimension > &gjit
          = geometry.jacobianInverseTransposed( local_point );

        typename LocalFunction::JacobianRangeType du, adu;
        uLocal.jacobian( quadrature[ quadraturePoint ], du );
        diffusion_operator_.diffusiveFlux( global_point, du, adu );

        for( unsigned int i = 0; i < numBaseFunctions; ++i )
        {
          typename BaseFunctionSet::JacobianRangeType j, dphi;
          baseSet.jacobian( i, quadrature[ quadraturePoint ], j );
          gjit.mv( j[ 0 ], dphi[ 0 ] );
          wLocal[ i ] += weight * (adu[ 0 ] * dphi[ 0 ]);
        }
      }
    }

    // boundary treatment
    const GridPart &gridPart = discreteFunctionSpace_.gridPart();
    for( Iterator it = discreteFunctionSpace_.begin(); it != end; ++it )
    {
      const Entity &entity = *it;
      if( !entity.hasBoundaryIntersections() )
        continue;

      const LagrangePointSet &lagrangePointSet = discreteFunctionSpace_.lagrangePointSet( entity );
      const LocalFunction uLocal = u.localFunction( entity );
      LocalFunction wLocal = w.localFunction( entity );

      const IntersectionIterator iend = gridPart.iend( entity );
      for( IntersectionIterator iit = gridPart.ibegin( entity ); iit != iend; ++iit )
      {
        const Intersection &intersection = *iit;
        if( !intersection.boundary() )
          continue;

        const int face = intersection.indexInInside();
        const FaceDofIterator fdend = lagrangePointSet.template endSubEntity< 1 >( face );
        for( FaceDofIterator fdit = lagrangePointSet.template beginSubEntity< 1 >( face ); fdit != fdend; ++fdit )
          wLocal[ *fdit ] = uLocal[ *fdit ];
      }
    }
#endif

  }


  template< class DiscreteFunctionImp, class DiffusionImp, class ReactionImp >
  template< class MatrixType >
  void DiscreteEllipticOperator< DiscreteFunctionImp, DiffusionImp, ReactionImp >::assemble_matrix ( MatrixType &global_matrix, bool boundary_treatment = true ) const
  {
    typedef typename MatrixType::LocalMatrixType LocalMatrix;

    global_matrix.reserve();
    global_matrix.clear();

    std::vector< typename BaseFunctionSet::JacobianRangeType > gradient_phi( discreteFunctionSpace_.mapper().maxNumDofs() );

    // micro scale base function:
    std::vector< RangeType > phi( discreteFunctionSpace_.mapper().maxNumDofs() );

    const Iterator end = discreteFunctionSpace_.end();
    for( Iterator it = discreteFunctionSpace_.begin(); it != end; ++it )
    {

      const Entity &entity = *it;
      const Geometry &geometry = entity.geometry();
      assert( entity.partitionType() == InteriorEntity );

      LocalMatrix local_matrix = global_matrix.localMatrix( entity, entity );

      const BaseFunctionSet &baseSet = local_matrix.domainBaseFunctionSet();
      const unsigned int numBaseFunctions = baseSet.numBaseFunctions();

      // for constant diffusion "2*discreteFunctionSpace_.order()" is sufficient, for the general case, it is better to use a higher order quadrature:
      Quadrature quadrature( entity, 2*discreteFunctionSpace_.order()+2 );
      const size_t numQuadraturePoints = quadrature.nop();
      for( size_t quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint )
      {
        // local (barycentric) coordinates (with respect to entity)
        const typename Quadrature::CoordinateType &local_point = quadrature.point( quadraturePoint );

        DomainType global_point = geometry.global( local_point );

        const double weight = quadrature.weight( quadraturePoint ) * geometry.integrationElement( local_point );

        // transposed of the the inverse jacobian
        const FieldMatrix< double, dimension, dimension > &inverse_jac
          = geometry.jacobianInverseTransposed( local_point );

        for( unsigned int i = 0; i < numBaseFunctions; ++i )
        {
          // jacobian of the base functions, with respect to the reference element
          typename BaseFunctionSet::JacobianRangeType gradient_phi_ref_element;
          baseSet.jacobian( i, quadrature[ quadraturePoint ], gradient_phi_ref_element );

          // multiply it with transpose of jacobian inverse to obtain the jacobian with respect to the real entity
          inverse_jac.mv( gradient_phi_ref_element[ 0 ], gradient_phi[ i ][ 0 ] );

          baseSet.evaluate( i, quadrature[ quadraturePoint ], phi[ i ]);

        }

        for( unsigned int i = 0; i < numBaseFunctions; ++i )
        {
          // A( \nabla \phi ) // diffusion operator evaluated in (x,\nabla \phi)
          typename LocalFunction::JacobianRangeType diffusion_in_gradient_phi;
          diffusion_operator_.diffusiveFlux( global_point, gradient_phi[ i ], diffusion_in_gradient_phi );
          for( unsigned int j = 0; j < numBaseFunctions; ++j )
           {
             local_matrix.add( j, i, weight * (diffusion_in_gradient_phi[ 0 ] * gradient_phi[ j ][ 0 ]) );

             if ( reaction_coefficient_ )
              {
                RangeType c;
                reaction_coefficient_->evaluate(global_point, c);
                local_matrix.add( j, i, weight * c * (phi[ i ][ 0 ] * phi[ j ][ 0 ]) );
              }

           }
        }
      }
    }


    // boundary treatment
    if ( boundary_treatment )
      {
        const GridPart &gridPart = discreteFunctionSpace_.gridPart();
        for( Iterator it = discreteFunctionSpace_.begin(); it != end; ++it )
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
  }

#if 1

  // assemble stiffness matrix for the jacobian matrix of the diffusion operator evaluated in the gradient of a certain discrete function (in case of the Newton method, it is the preceeding iterate u_H^{(n-1)} )
  // stiffness matrix with entries
  // \int JA(\nabla disc_func) \nabla phi_i \nabla phi_j 
  // (here, JA denotes the jacobian matrix of the diffusion operator A)
  template< class DiscreteFunctionImp, class DiffusionImp, class ReactionImp >
  template< class MatrixType >
  void DiscreteEllipticOperator< DiscreteFunctionImp, DiffusionImp, ReactionImp >::assemble_jacobian_matrix ( DiscreteFunction &disc_func, MatrixType &global_matrix, bool boundary_treatment = true  ) const
  {
    typedef typename MatrixType::LocalMatrixType LocalMatrix;

    typedef typename DiscreteFunction :: LocalFunctionType
      LocalFunction;

    global_matrix.reserve();
    global_matrix.clear();

    std::vector< typename BaseFunctionSet::JacobianRangeType > gradient_phi( discreteFunctionSpace_.mapper().maxNumDofs() );

    // micro scale base function:
    std::vector< RangeType > phi( discreteFunctionSpace_.mapper().maxNumDofs() );

    const Iterator end = discreteFunctionSpace_.end();
    for( Iterator it = discreteFunctionSpace_.begin(); it != end; ++it )
    {

      const Entity &entity = *it;
      const Geometry &geometry = entity.geometry();
      assert( entity.partitionType() == InteriorEntity );

      LocalMatrix local_matrix = global_matrix.localMatrix( entity, entity );
      LocalFunction local_disc_function = disc_func.localFunction( entity ); 

      const BaseFunctionSet &baseSet = local_matrix.domainBaseFunctionSet();
      const unsigned int numBaseFunctions = baseSet.numBaseFunctions();

      // for constant diffusion "2*discreteFunctionSpace_.order()" is sufficient, for the general case, it is better to use a higher order quadrature:
      Quadrature quadrature( entity, 2*discreteFunctionSpace_.order()+2 );
      const size_t numQuadraturePoints = quadrature.nop();
      for( size_t quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint )
      {
        // local (barycentric) coordinates (with respect to entity)
        const typename Quadrature::CoordinateType &local_point = quadrature.point( quadraturePoint );

        DomainType global_point = geometry.global( local_point );

        const double weight = quadrature.weight( quadraturePoint ) * geometry.integrationElement( local_point );

        // transposed of the the inverse jacobian
        const FieldMatrix< double, dimension, dimension > &inverse_jac
          = geometry.jacobianInverseTransposed( local_point );

        for( unsigned int i = 0; i < numBaseFunctions; ++i )
        {
          // jacobian of the base functions, with respect to the reference element
          typename BaseFunctionSet::JacobianRangeType gradient_phi_ref_element;
          baseSet.jacobian( i, quadrature[ quadraturePoint ], gradient_phi_ref_element );

          // multiply it with transpose of jacobian inverse to obtain the jacobian with respect to the real entity
          inverse_jac.mv( gradient_phi_ref_element[ 0 ], gradient_phi[ i ][ 0 ] );

          baseSet.evaluate( i, quadrature[ quadraturePoint ], phi[ i ]);

        }

        for( unsigned int i = 0; i < numBaseFunctions; ++i )
        {

          typename BaseFunctionSet::JacobianRangeType grad_local_disc_function;
          local_disc_function.jacobian(quadrature[ quadraturePoint ], grad_local_disc_function);
          // here: no multiplication with jacobian inverse transposed required!

          // JA( \nabla u_H ) \nabla phi_i // jacobian of diffusion operator evaluated in (x,grad_local_disc_function) in direction of the gradient of the current base function
          typename LocalFunction::JacobianRangeType jac_diffusion_flux;
          diffusion_operator_.jacobianDiffusiveFlux( global_point, grad_local_disc_function, gradient_phi[ i ], jac_diffusion_flux );

          for( unsigned int j = 0; j < numBaseFunctions; ++j )
            {
             local_matrix.add( j, i, weight * (jac_diffusion_flux[ 0 ] * gradient_phi[ j ][ 0 ]) );

             if ( reaction_coefficient_ )
              {
                RangeType c;
                reaction_coefficient_->evaluate(global_point, c);
                local_matrix.add( j, i, weight * c * (phi[ i ][ 0 ] * phi[ j ][ 0 ]) );
              }
            }
        }
      }
    }

    // boundary treatment
    if ( boundary_treatment )
      {
        const GridPart &gridPart = discreteFunctionSpace_.gridPart();
        for( Iterator it = discreteFunctionSpace_.begin(); it != end; ++it )
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


  }


#endif

//!NOTE DEPRECATED !!!!!!!!!!!!!!!!!!!!!!!! (loeschen):
#if 0
  // assembleRHS
  // -----------

  template< class Function, class DiscreteFunction >
  void assembleRHS ( const Function &function, DiscreteFunction &rhs )
  {
    typedef typename DiscreteFunction::DiscreteFunctionSpaceType DiscreteFunctionSpace;
    typedef typename DiscreteFunction::LocalFunctionType LocalFunction;

    typedef typename DiscreteFunctionSpace::BaseFunctionSetType BaseFunctionSet;
    typedef typename DiscreteFunctionSpace::IteratorType Iterator;
    typedef typename Iterator::Entity Entity;
    typedef typename Entity::Geometry Geometry;

    typedef typename DiscreteFunctionSpace::GridPartType GridPart;
    typedef CachingQuadrature< GridPart, 0 > Quadrature;

    const DiscreteFunctionSpace &discreteFunctionSpace = rhs.space();

    const Iterator end = discreteFunctionSpace.end();
    for( Iterator it = discreteFunctionSpace.begin(); it != end; ++it )
    {
      const Entity &entity = *it;
      const Geometry &geometry = entity.geometry();
      assert( entity.partitionType() == InteriorEntity );

      LocalFunction rhsLocal = rhs.localFunction( entity );
      
      const BaseFunctionSet &baseSet = rhsLocal.baseFunctionSet();
      const unsigned int numBaseFunctions = baseSet.numBaseFunctions();
           
      Quadrature quadrature( entity, 2*discreteFunctionSpace.order()+1 );
      const size_t numQuadraturePoints = quadrature.nop();
      for( size_t quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint )
      {
        const typename Quadrature::CoordinateType &x = quadrature.point( quadraturePoint );
        const double weight = quadrature.weight( quadraturePoint ) * geometry.integrationElement( x );

        typename Function::RangeType f;
        function.evaluate( geometry.global( x ), f );

        for( unsigned int i = 0; i < numBaseFunctions; ++i )
        {
          typename BaseFunctionSet::RangeType phi;
          baseSet.evaluate( i, quadrature[ quadraturePoint ], phi );
          rhsLocal[ i ] += weight * (f * phi);
        }
      }
    }
  }



  // boundaryTreatment
  // -----------------

  template< class GridFunction, class DiscreteFunction >
  void boundaryTreatment ( const GridFunction &exactSolution,
                           DiscreteFunction &rhs, DiscreteFunction &solution )
  {
    typedef typename DiscreteFunction::DiscreteFunctionSpaceType DiscreteFunctionSpace;
    typedef typename DiscreteFunction::LocalFunctionType LocalFunction;

    typedef typename GridFunction::LocalFunctionType LocalExactSolution;

    typedef typename DiscreteFunctionSpace::IteratorType Iterator;
    typedef typename Iterator::Entity Entity;

    typedef typename DiscreteFunctionSpace::GridPartType GridPart;
    typedef typename GridPart::IntersectionIteratorType IntersectionIterator;
    typedef typename IntersectionIterator::Intersection Intersection;

    typedef typename DiscreteFunctionSpace::LagrangePointSetType LagrangePointSet;
    typedef typename LagrangePointSet::template Codim< 1 >::SubEntityIteratorType FaceDofIterator;

    const DiscreteFunctionSpace &discreteFunctionSpace = rhs.space();

    const Iterator end = discreteFunctionSpace.end();
    for( Iterator it = discreteFunctionSpace.begin(); it != end; ++it )
    {
      const Entity &entity = *it;
      if( !entity.hasBoundaryIntersections() )
        continue;

      LocalExactSolution exactLocal = exactSolution.localFunction( entity );
      LocalFunction rhsLocal = rhs.localFunction( entity );
      LocalFunction solutionLocal = solution.localFunction( entity );

      const LagrangePointSet &lagrangePointSet = discreteFunctionSpace.lagrangePointSet( entity );

      const IntersectionIterator iend = discreteFunctionSpace.gridPart().iend( entity );
      for( IntersectionIterator iit = discreteFunctionSpace.gridPart().ibegin( entity ); iit != iend; ++iit )
      {
        const Intersection &intersection = *iit;
        if( !intersection.boundary() )
          continue;

        const int face = intersection.indexInInside();
        const FaceDofIterator fdend = lagrangePointSet.template endSubEntity< 1 >( face );
        for( FaceDofIterator fdit = lagrangePointSet.template beginSubEntity< 1 >( face ); fdit != fdend; ++fdit )
        {
          typename LocalExactSolution::RangeType phi;
          exactLocal.evaluate( lagrangePointSet[ *fdit ], phi );
          rhsLocal[ *fdit ] = phi[ 0 ];
          solutionLocal[ *fdit ] = phi[ 0 ];
        }
      }
    }
  }
#endif

}

#endif // #ifndef DiscreteElliptic_HH
