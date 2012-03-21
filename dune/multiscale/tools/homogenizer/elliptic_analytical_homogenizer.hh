
//das alles klappt (mathematisch) nur in 2-D!!! Für Tensoren, die:
// 1. irgendwelche Elliptizitätsbedingungen erfüllen
// 2. A(x,y) = A(y)
// 3. a_i_j(y) = a_i_j(y_1,y_2) = a_i_j(y_1)
// 4. symmetrisch

#ifndef DUNE_ANALYTICALHOMOGENIZER_HH
#define DUNE_ANALYTICALHOMOGENIZER_HH

// where the quadratures are defined 
#include <dune/fem/quadrature/cachingquadrature.hh>

namespace Dune 
{
  
  template < class GridImp, class TensorImp > 
  class AnalyticalHomogenizer
  {
   typedef GridImp GridType;

   typedef LeafGridPart< GridType > GridPartType;

   typedef FunctionSpace < double , double , 2 , 1 > FunctionSpaceType; 

   typedef LagrangeDiscreteFunctionSpace < FunctionSpaceType, GridPartType, 1 >
      DiscreteFunctionSpaceType;

   typedef AdaptiveDiscreteFunction < DiscreteFunctionSpaceType > DiscreteFunctionType;

   typedef typename FunctionSpaceType::DomainType DomainType; 

   typedef typename FunctionSpaceType::RangeType RangeType;

    typedef TensorImp TensorType;

    typedef typename DiscreteFunctionSpaceType::IteratorType
       IteratorType;

    typedef typename GridType::template Codim<0>::Entity
       EntityType; 
    
    typedef typename GridType::template Codim<0>::Geometry
       EnGeometryType; 
    
    typedef typename EntityType::ctype 
       coordType; 
      
    enum { dimension = GridType::dimension};
    enum { spacePolOrd = DiscreteFunctionSpaceType :: polynomialOrder };

    mutable FieldMatrix< RangeType, dimension, dimension > a;


    // dgf file that describes the perforated domain
    std :: string &filename_;

  public:

    AnalyticalHomogenizer( std :: string &filename )
    : filename_( filename )
    {
    }

    FieldMatrix< RangeType, dimension, dimension > getHomTensor
         ( const TensorType &tensor, int polOrd = (2 * spacePolOrd + 2) ) const
    {

      std :: cout << "WARNING! Use of deprecated homogenizer, which requires deprecated use of 'evaluate' in Diffusion-class!" << std :: endl; 

      GridPtr< GridType > gridptr( filename_ ); 

      gridptr->globalRefine( 12 );

      GridPartType gridPart( *gridptr );

      DiscreteFunctionSpaceType discreteFunctionSpace( gridPart );

      FieldMatrix< RangeType, dimension, dimension > tensorHom(0.0);

      IteratorType endit = discreteFunctionSpace.end();
      for(IteratorType it = discreteFunctionSpace.begin(); it != endit ; ++it)
       {
        // entity
        const EntityType& entity = *it;

        // create quadrature for given geometry type 
        CachingQuadrature <GridPartType , 0 > quadrature(entity,polOrd); 

        // get geoemetry of entity
        const EnGeometryType& geometry = entity.geometry();

        // integrate 
        const int quadratureNop = quadrature.nop();
        for(int quadraturePoint = 0; quadraturePoint < quadratureNop; ++quadraturePoint)
         {
          const double det = quadrature.weight(quadraturePoint) * 
          geometry.integrationElement(quadrature.point(quadraturePoint));

          tensor.evaluate( 0, 0,
                           geometry.global( quadrature.point( quadraturePoint ) ),
                           a[ 0 ][ 0 ] );

          tensorHom[ 0 ][ 0 ] += det * ( 1 / a[ 0 ][ 0 ] );
         }
       } 

      tensorHom[ 0 ][ 0 ] = 1 / tensorHom[ 0 ][ 0 ];

      for(IteratorType it = discreteFunctionSpace.begin(); it != endit ; ++it)
       {
        // entity
        const EntityType& entity = *it;

        // create quadrature for given geometry type 
        CachingQuadrature <GridPartType , 0 > quadrature(entity,polOrd); 

        // get geoemetry of entity
        const EnGeometryType& geometry = entity.geometry();

        // integrate 
        const int quadratureNop = quadrature.nop();
        for(int quadraturePoint = 0; quadraturePoint < quadratureNop; ++quadraturePoint)
         {
          const double det = quadrature.weight(quadraturePoint) * 
              geometry.integrationElement(quadrature.point(quadraturePoint));

          tensor.evaluate( 0, 0,
                           geometry.global( quadrature.point( quadraturePoint ) ),
                           a[ 0 ][ 0 ] );

          tensor.evaluate( 0, 1,
                           geometry.global( quadrature.point( quadraturePoint ) ),
                           a[ 0 ][ 1 ] );

          tensorHom[ 0 ][ 1 ] += det * ( a[ 0 ][ 1 ] / a[ 0 ][ 0 ] );
         }
       }

       tensorHom[ 0 ][ 1 ] *= tensorHom[ 0 ][ 0 ];

       tensorHom[ 1 ][ 0 ] = tensorHom[ 0 ][ 1 ];


       for(IteratorType it = discreteFunctionSpace.begin(); it != endit ; ++it)
        {
          // entity
          const EntityType& entity = *it;

          // create quadrature for given geometry type 
          CachingQuadrature <GridPartType , 0 > quadrature(entity,polOrd); 

          // get geoemetry of entity
          const EnGeometryType& geometry = entity.geometry();

          // integrate 
          const int quadratureNop = quadrature.nop();
          for(int quadraturePoint = 0; quadraturePoint < quadratureNop; ++quadraturePoint)
           {
             const double det = quadrature.weight(quadraturePoint) * 
                   geometry.integrationElement(quadrature.point(quadraturePoint));

             tensor.evaluate( 0, 0,
                              geometry.global( quadrature.point( quadraturePoint ) ),
                              a[ 0 ][ 0 ] );

             tensor.evaluate( 0, 1,
                              geometry.global( quadrature.point( quadraturePoint ) ),
                              a[ 0 ][ 1 ] );

             tensor.evaluate( 1, 0,
                              geometry.global( quadrature.point( quadraturePoint ) ),
                              a[ 1 ][ 0 ] );

             tensor.evaluate( 1, 1,
                              geometry.global( quadrature.point( quadraturePoint ) ),
                              a[ 1 ][ 1 ] );

             tensorHom[ 1 ][ 1 ] += det * 
                                     ( a[ 1 ][ 1 ] - ( (a[ 0 ][ 1 ] * a[ 1 ][ 0 ]) / a[ 0 ][ 0 ] ) );
            }
         }

       tensorHom[ 1 ][ 1 ] += (tensorHom[ 1 ][ 0 ] * tensorHom[ 0 ][ 1 ]) /
                                    tensorHom[ 0 ][ 0 ];

       std :: cout << "analytical: A_homogenized[0][0] = " << tensorHom[0][0] << std :: endl;
       std :: cout << "analytical: A_homogenized[0][1] = " << tensorHom[0][1] << std :: endl;
       std :: cout << "analytical: A_homogenized[1][0] = " << tensorHom[1][0] << std :: endl;
       std :: cout << "analytical: A_homogenized[1][1] = " << tensorHom[1][1] << std :: endl;

       return tensorHom;

    } // end of method

}; // end of class



} // end namespace 
#endif
