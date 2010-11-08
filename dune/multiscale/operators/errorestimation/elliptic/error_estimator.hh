#ifndef DUNE_ERRORESTIMATER_HH
#define DUNE_ERRORESTIMATER_HH

// where the quadratures are defined 
#include <dune/fem/quadrature/cachingquadrature.hh>
#include "../../cell_problem_solving/elliptic/old_cellproblemsolver.hh"

//! ########### ELLIPTIC ###################################################

namespace Dune 
{

//! Klasse loeschen?:


  template< class DiscreteFunctionImp >
  class DiscFuncAdapter
  {
  public:
    
    typedef DiscreteFunctionImp DiscreteFunctionType;
  
    typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType
       DiscreteFunctionSpaceType;
    
    typedef typename DiscreteFunctionType::RangeType
       RangeType;
    
    typedef typename DiscreteFunctionSpaceType::IteratorType
       IteratorType;
    
    typedef typename DiscreteFunctionSpaceType::GridType
       GridType;
    
    typedef typename DiscreteFunctionSpaceType::GridPartType
       GridPartType;
    
    typedef typename GridType::template Codim<0>::Entity
       EntityType; 
    
    typedef typename GridType::template Codim<0>::Geometry
       EnGeometryType; 
    
    typedef typename EntityType::ctype 
       coordType; 
    
    typedef typename DiscreteFunctionType::LocalFunctionType
       LocalFunctionType;

    enum { dimension = GridType::dimension};
    enum { spacePolOrd = DiscreteFunctionSpaceType :: polynomialOrder }; 

  public:
   
    //discreteFunction = cof * discreteFunction
    //template <class FunctionType> 
    static void multScalar(DiscreteFunctionType &discreteFunction, double &scalar)
    {
      typedef typename DiscreteFunctionType :: DofIteratorType DofIteratorType;

      const DofIteratorType end = discreteFunction.dend();
      for( DofIteratorType it = discreteFunction.dbegin(); it != end; ++it )
        *it = (*it) * scalar;
    }

    //discreteFunction1 = discreteFunction1 + discreteFunction2
    static void add(DiscreteFunctionType &discreteFunction1,
                    DiscreteFunctionType &discreteFunction2)
    {

     typedef typename DiscreteFunctionType :: DofIteratorType DofIteratorType;


     double value[10000];
     for( int k = 0; k < 10000; ++k )
       {
        value[ k ] = 0;
       }

     int globalDofNumber = 0;


     const DofIteratorType end2 = discreteFunction2.dend();
     for( DofIteratorType it = discreteFunction2.dbegin(); it != end2; ++it )
       {
        value[ globalDofNumber ] = (*it);
        globalDofNumber += 1;
       }


      globalDofNumber = 0;
      const DofIteratorType end1 = discreteFunction1.dend();
      for( DofIteratorType it = discreteFunction1.dbegin(); it != end1; ++it )
        {
          *it = value[ globalDofNumber ] + (*it);
          globalDofNumber += 1;
        }
    }
 
    static void setZero(DiscreteFunctionType &discreteFunction)
    {
      typedef typename DiscreteFunctionType :: DofIteratorType DofIteratorType;

      const DofIteratorType end = discreteFunction.dend();
      for( DofIteratorType it = discreteFunction.dbegin(); it != end; ++it )
        {
         *it = 0;
        }
    }

  };





  template <class PeriodicDiscreteFunctionImp,
            class DiscreteFunctionImp,
            class TensorImp,
            class MassTermImp,
            class RHSFunctionImp > 
//NOTE: das zweite Argument DiscreteFunctionImp sobald wie moeglich wieder loeschen. Alle 'Ableitungen' und Variablen von DiscreteFunctionImp ebenfalls loeschen und die zugehörigen von PeriodicDiscreteFunctionImp ersetzen. Aktuell ist es nur da, weil verschiedene Dinge wie der IntersectionIterator für das periodic Gridpart noch nicht implementiert sind.
  class ErrorEstimator
  {

    typedef TensorImp TensorType;
    typedef MassTermImp MassType;
    typedef RHSFunctionImp RHSFunctionType;

    //! Necessary typedefs for the PeriodicDiscreteFunctionImp:

    typedef PeriodicDiscreteFunctionImp PeriodicDiscreteFunctionType;

    typedef typename PeriodicDiscreteFunctionType       :: LocalFunctionType PeriodicLocalFunctionType;
    typedef typename PeriodicDiscreteFunctionType       :: FunctionSpaceType PeriodicDiscreteFunctionSpaceType;
    typedef typename PeriodicDiscreteFunctionSpaceType  :: GridPartType      PeriodicGridPartType;
    typedef typename PeriodicDiscreteFunctionSpaceType  :: GridType          PeriodicGridType;
    typedef typename PeriodicDiscreteFunctionType       :: RangeType         RangeType;
    typedef typename PeriodicDiscreteFunctionType       :: RangeFieldType    TimeType;
    typedef typename PeriodicDiscreteFunctionType       :: DomainType        DomainType;
    typedef typename PeriodicDiscreteFunctionSpaceType  :: JacobianRangeType PeriodicJacobianRangeType;
    typedef typename PeriodicDiscreteFunctionSpaceType  :: IteratorType      PeriodicIteratorType;

    typedef typename PeriodicGridPartType                  :: IntersectionIteratorType PeriodicIntersectionIteratorType ;
    typedef typename PeriodicGridType :: template Codim<0> :: Entity                   PeriodicEntityType; 
    typedef typename PeriodicGridType :: template Codim<0> :: EntityPointer            PeriodicEntityPointerType; 
    typedef typename PeriodicGridType :: template Codim<0> :: Geometry                 PeriodicEntityGeometryType; 
    typedef typename PeriodicGridType :: template Codim<1> :: Geometry                 PeriodicFaceGeometryType;

    typedef CachingQuadrature < PeriodicGridPartType , 0 > PeriodicEntityQuadratureType;
    typedef CachingQuadrature < PeriodicGridPartType , 1 > PeriodicFaceQuadratureType;
    //! Necessary typedefs for the DiscreteFunctionImp:

    typedef DiscreteFunctionImp DiscreteFunctionType;

    typedef typename DiscreteFunctionType      :: LocalFunctionType LocalFunctionType;
    typedef typename DiscreteFunctionType      :: FunctionSpaceType DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctionType      :: DofIteratorType   DofIteratorType;
    typedef typename DiscreteFunctionSpaceType :: GridPartType      GridPartType;
    typedef typename DiscreteFunctionSpaceType :: GridType          GridType;
    typedef typename DiscreteFunctionSpaceType :: JacobianRangeType JacobianRangeType;
    typedef typename DiscreteFunctionSpaceType :: IteratorType      IteratorType;

    typedef typename GridPartType                  :: IntersectionIteratorType IntersectionIteratorType;
    typedef typename GridType :: template Codim<0> :: Entity                   EntityType; 
    typedef typename GridType :: template Codim<0> :: EntityPointer            EntityPointerType; 
    typedef typename GridType :: template Codim<0> :: Geometry                 EntityGeometryType; 
    typedef typename GridType :: template Codim<1> :: Geometry                 FaceGeometryType;
    typedef typename DiscreteFunctionSpaceType     :: BaseFunctionSetType      BaseFunctionSetType;

    typedef CachingQuadrature < GridPartType , 0 > EntityQuadratureType;
    typedef CachingQuadrature < GridPartType , 1 > FaceQuadratureType;

    enum { dimension = PeriodicGridType :: dimension};
    enum { spacePolOrd = PeriodicDiscreteFunctionSpaceType :: polynomialOrder }; 
    enum { maxnumOfBaseFct = 100 }; 

  private:

    const PeriodicDiscreteFunctionSpaceType &periodicDiscreteFunctionSpace_;
    const DiscreteFunctionSpaceType &discreteFunctionSpace_;
    const DiscreteFunctionSpaceType &auxiliaryDiscreteFunctionSpace_;
    // an auxiliaryDiscreteFunctionSpace to get an Intersection Iterator for the periodicDiscreteFunctionSpace
    // (for the periodic grid partition there is no usable intersection iterator implemented, therefore we use the intersection iterator for the corresponding non-periodic grid partition (this not efficient and increases the estimated error, but it works)


    const TensorType &tensor_;
    const MassType &massTerm_;
    const RHSFunctionType &rhs_;

    const RangeType epsilon_est_;
    std :: string &cell_solution_location_;

  public:

    ErrorEstimator ( const PeriodicDiscreteFunctionSpaceType &periodicDiscreteFunctionSpace,
                     const DiscreteFunctionSpaceType &discreteFunctionSpace,
                     const DiscreteFunctionSpaceType &auxiliaryDiscreteFunctionSpace,
                     const TensorType &tensor,
                     const MassType &massTerm,
                     const RHSFunctionType &rhs,
                     const RangeType &epsilon_est,
                           std :: string &cell_solution_location )
    : periodicDiscreteFunctionSpace_( periodicDiscreteFunctionSpace ),
      discreteFunctionSpace_( discreteFunctionSpace ),
      auxiliaryDiscreteFunctionSpace_( auxiliaryDiscreteFunctionSpace ),
      tensor_( tensor ),
      massTerm_( massTerm ),
      rhs_ ( rhs ),
      epsilon_est_( epsilon_est ),
      cell_solution_location_( cell_solution_location )
    {
    }

   int get_cell_problem_id( const EntityPointerType &it_now, int &base_func_id )
   {

    int global_key = 0;

    int counter = 0;

    IteratorType endit = discreteFunctionSpace_.end();

    for(IteratorType it = discreteFunctionSpace_.begin(); it != endit ; ++it)
     {

        // entity
        const EntityType& entity = *it;

        const BaseFunctionSetType baseSet
            = discreteFunctionSpace_.baseFunctionSet( entity );

        // number of base functions on entity
        const int numBaseFunctions = baseSet.numBaseFunctions();

        for( int i = 0; i < numBaseFunctions; ++i )
          {
            if ( (it == it_now) && ( i==base_func_id ) )
             {global_key = counter;}

            counter = counter + 1;
          }

     }
    return global_key;
   }

// the old method:
#if 0
    // method that gets the local mesh size H_entity
    RangeType getH( const EntityType &entity)
    {
      const GridPartType &gridPart = discreteFunctionSpace_.gridPart();

      EntityQuadratureType entityQuadrature(entity,0); //PolOrd = 0

      const EntityGeometryType& geometry = entity.geometry();
 
      const RangeType entityVolume = entityQuadrature.weight(0) * 
              geometry.integrationElement(entityQuadrature.point(0));

      RangeType faceVolume(0.0);
      int numberOfFaces = 0;

      // get the 'volume' of every face of the entity and devide the result by the number of faces to get an averaged 'faceVolume' 
      IntersectionIteratorType endnit = gridPart.iend(entity);
      for(IntersectionIteratorType nit = gridPart.ibegin(entity); nit != endnit ; ++nit)
        {
          FaceQuadratureType innerFaceQuadrature(gridPart, nit, 0 , FaceQuadratureType::INSIDE); 

          DomainType scaledOuterNormal =
               nit.integrationOuterNormal(innerFaceQuadrature.localPoint(0));

          // get 'volume' of the visited face (this only works because we do not have curved faces):
          RangeType visitedFaceVolume(0.0);
          for ( int k = 0; k < dimension; ++k ) 
                visitedFaceVolume += scaledOuterNormal[k] * scaledOuterNormal[k];
          visitedFaceVolume = sqrt(visitedFaceVolume);

          faceVolume += visitedFaceVolume;

          numberOfFaces += 1;
        }

       faceVolume = faceVolume / numberOfFaces;


       // instead of using h_entity = entityVolume, we set:
       RangeType h_entity = entityVolume / faceVolume;

       return h_entity;

    }
#endif
// the new method:
    //! method to get the local mesh size H_entity (of the macro mesh)
    // works only for our 2D examples!!!! 
    RangeType getH( const EntityType &entity)
    {

      // entity_H means H (the diameter of the entity)
      RangeType entity_H = 0.0;

      const GridPartType &gridPart = discreteFunctionSpace_.gridPart();
 
      // compute the size of the faces of the entities and selected the largest.
      IntersectionIteratorType endnit = gridPart.iend(entity);
      for(IntersectionIteratorType nit = gridPart.ibegin(entity); nit != endnit ; ++nit)
        {
          FaceQuadratureType innerFaceQuadrature(gridPart, nit, 0 , FaceQuadratureType::INSIDE); 

          DomainType scaledOuterNormal =
               nit.integrationOuterNormal(innerFaceQuadrature.localPoint(0));

          // get 'volume' of the visited face (this only works because we do not have curved faces):
          RangeType visitedFaceVolume(0.0);
          for ( int k = 0; k < dimension; ++k ) 
                visitedFaceVolume += scaledOuterNormal[k] * scaledOuterNormal[k];
          visitedFaceVolume = sqrt(visitedFaceVolume);

          if (visitedFaceVolume > entity_H)
            entity_H = visitedFaceVolume;

        }

       return entity_H;

    }

    // return ||f||_{L^2(T)}
    RangeType errorf( const EntityType &entity )
    {

        // create quadrature for given geometry type 
        CachingQuadrature <GridPartType , 0 > entityQuadrature(entity, spacePolOrd ); 

        // get geoemetry of entity
        const EntityGeometryType& geometry = entity.geometry();

        RangeType y(0);
        RangeType localError(0);

        const int quadratureNop = entityQuadrature.nop();
        for(int quadraturePoint = 0; quadraturePoint < quadratureNop; ++quadraturePoint)
        {
          const double det = entityQuadrature.weight(quadraturePoint) * 
              geometry.integrationElement(entityQuadrature.point(quadraturePoint));

          rhs_.evaluate(geometry.global( entityQuadrature.point(quadraturePoint) ),y);
          y = y * y;

          localError += det * y;
        }

        localError = sqrt(localError);

        return localError;
    }
 

//calculate:
// \sum_(all faces of Y) 
//   \int_E h³_entity ([( A_h(x_j,y) (\nabla_x u_H(x_j) + \nabla_y (K_h(u_H))]_E * n)²
    RangeType errorXi
             (const DomainType &globalPoint,
              PeriodicDiscreteFunctionType &corrector_u_H,
              JacobianRangeType &grad_u_H,
              int polOrd = 0 /*(2 * spacePolOrd + 2)*/) const
    {

    //! WICHTIGE NOTIZ: das hier ist natürlich alles suboptimal, da wir statt dem grid auf dem periodic disrecte function space das grid auf dem normalen discrete function space verwenden. leider lässt sich das momentan noch nicht umgehen. Damit ergeben sich zwei Problem:
    //!  1. H = h (da wir das MacroGrid hier als MicroGrid verwenden --> ist behoben
    //!  2. wir verschenken die Rand Sprünge

     bool homogeneous_structure = tensor_.homogeneous();

     const GridPartType &gridPart = auxiliaryDiscreteFunctionSpace_.gridPart();

     RangeType result(0.0);

     IteratorType endit = auxiliaryDiscreteFunctionSpace_.end();
     for(IteratorType it = auxiliaryDiscreteFunctionSpace_.begin(); it != endit ; ++it)
      { //start entity iterator

// das folgende nutzt merfach direkt, dass bei LagrangeSpace mit PolOrder 1, die Gradienten auf jeder Entity konstant sind und es deshalb egal ist, welcher Quadraturpunkt verwendet wird. Wir verwenden für EntityQuadraturen stets die Mittelpunkt Quadratur. Damit ist der einzige Quadraturpunkt stets der Mittelpunkt der Entity: entityQuadrature, 0 (0 = quadPoint).
// Entsprechend gibt es auch keine Schleife der Art:

//    const int quadratureNop = entityQuadrature.nop();
//    for(int entityQuadPoint = 0; entityQuadPoint < quadratureNop; ++entityQuadPoint)

       // entity
       const EntityType& entity = *it;
     
       // create quadrature for given geometry type 
       EntityQuadratureType entityQuadrature(entity,0); //PolOrd = 0  
       // we also need a quadrature for all the neighbor entities,
       // they are implemented in the next loop

       // get geoemetry of the entity
       const EntityGeometryType& geometry = entity.geometry();
       // we also need the geometry of all the neighbor entities,
       // they are implemented in the next loop
 
       const RangeType entityVolume = entityQuadrature.weight(0) * 
              geometry.integrationElement(entityQuadrature.point(0));

       IntersectionIteratorType endnit = gridPart.iend(entity);
       for(IntersectionIteratorType nit = gridPart.ibegin(entity); nit != endnit ; ++nit)
          {
           //the following distinction is redundant if we are able to use the intersection iterator and face quadratures for the periodic gridpart:

           // if the face is not part of the outer boundary
           if ( nit.neighbor() ) //if there is a neighbor entity
           {

           //to save gradient of Corrector(Phi_H) on entity:
           JacobianRangeType gradCor(0.0);
           //to save gradient of Corrector(Phi_H) on neighbor entity:
           JacobianRangeType gradCorNeigh(0.0);

            //?const FaceGeometryType& faceGeometry = *nit->intersectionGlobal(); //geometry();

           // get the neighbor entity:
           const EntityPointerType neighborEntityPointer = nit.outside();
           const EntityType& neighborEntity = *neighborEntityPointer;
           // get quadrature of neighbor entity
           EntityQuadratureType neighborEntityQuadrature(neighborEntity,0); //PolOrd = 0 
           // only one qudraturePoint!
           const EntityGeometryType& neighborGeometry = neighborEntity.geometry();

           const PeriodicLocalFunctionType neighborLocalFunction = corrector_u_H.localFunction(neighborEntity);
           const PeriodicLocalFunctionType localFunction = corrector_u_H.localFunction(entity);  

           RangeType a_inside[ dimension ][ dimension ];
           if (homogeneous_structure == true )
            {
             for( int k = 0; k < dimension; ++k )
               for( int l = 0; l < dimension; ++l ) 
               {
                tensor_.evaluate( l, k,
                               globalPoint,
                               geometry.global( entityQuadrature.point( 0 ) ),
                               a_inside[ l ][ k ] );
               }
             }
           else
            {
             DomainType y_eps; // x_j + \eps_{est} * y
             for( int k = 0; k < dimension; ++k )
                y_eps[ k ] = globalPoint[ k ] + epsilon_est_ * geometry.global( entityQuadrature.point( 0 ) )[ k ];

             for( int k = 0; k < dimension; ++k )
               for( int l = 0; l < dimension; ++l ) 
               {
                tensor_.evaluate( l, k, y_eps, a_inside[ l ][ k ] );
               }
             }

           //to save A_h(EntityPoint) \nabla Phi_H(globalPoint):
           RangeType w1_inside[dimension];
           for( int k = 0; k < dimension; ++k ) 
               w1_inside[ k ] = 0; 

           for( int k = 0; k < dimension; ++k )
             for( int l = 0; l < dimension; ++l )
                w1_inside[ k ] += a_inside[ k ][ l ] * grad_u_H[ 0 ][ l ];


           RangeType a_outside[ dimension ][ dimension ];
           if (homogeneous_structure == true )
            {
             for( int k = 0; k < dimension; ++k )
               for( int l = 0; l < dimension; ++l ) 
                {
                 tensor_.evaluate( l, k,
                                  globalPoint,
                                  neighborGeometry.global( neighborEntityQuadrature.point( 0 ) ),
                                  a_outside[ l ][ k ] );
                }
             }
            else
             {
              DomainType y_eps; // x_j + \eps_{est} * y
              for( int k = 0; k < dimension; ++k )
                 y_eps[ k ] = globalPoint[ k ] + epsilon_est_ * neighborGeometry.global( neighborEntityQuadrature.point( 0 ) )[ k ];

              for( int k = 0; k < dimension; ++k )
                for( int l = 0; l < dimension; ++l ) 
                 {
                  tensor_.evaluate( l, k, y_eps, a_outside[ l ][ k ] );
                 }
              }

            //to save A_h(neighborEntityPoint) \nabla Phi_H(globalPoint):
            RangeType w1_outside[dimension];
            for( int k = 0; k < dimension; ++k ) 
                w1_outside[ k ] = 0; 

            for( int k = 0; k < dimension; ++k )
              for( int l = 0; l < dimension; ++l )
                 w1_outside[ k ] += a_outside[ k ][ l ] * grad_u_H[ 0 ][ l ];

           FaceQuadratureType innerFaceQuadrature(gridPart, nit, polOrd, FaceQuadratureType::INSIDE);
           // we could also use an outerFaceQuadrature:
           // FaceQuadratureType outerFaceQuadrature(gridPart, nit, polOrd, FaceQuadratureType::OUTSIDE);
           // Note that the results would be the same

           // innerFaceQuadrature and outerFaceQuadrature are equivalent, therefor we it is sufficient to stay restricted to the innerFaceQuadrature
           const int quadratureNop = innerFaceQuadrature.nop();
           for(int faceQuadPoint = 0; faceQuadPoint < quadratureNop; ++faceQuadPoint)
             {

              DomainType unitOuterNormal =
                nit.unitOuterNormal(innerFaceQuadrature.localPoint(faceQuadPoint));

              //outer normal Scaled with face size:
              DomainType scaledOuterNormal =
                nit.integrationOuterNormal(innerFaceQuadrature.localPoint(faceQuadPoint));

              // get 'volume' of the face (this only works because we do not have curved faces):
              RangeType faceVolume(0.0);
              for ( int k = 0; k < dimension; ++k ) 
                faceVolume += scaledOuterNormal[k] * scaledOuterNormal[k];
              faceVolume = sqrt(faceVolume);

              // instead of using h_entity = entityVolume, we set:
              RangeType h_entity = entityVolume / faceVolume;
              h_entity = faceVolume; //! loeschen?
              // then the implementation is independent of the dimension

              localFunction.jacobian( innerFaceQuadrature , faceQuadPoint , gradCor );
              neighborLocalFunction.jacobian( innerFaceQuadrature , faceQuadPoint , gradCorNeigh );

              //to save A_h(EntityPoint) \nabla Corrector(Phi_H)(EntityQuadPoint).
              RangeType w2_inside[dimension];
              for( int k = 0; k < dimension; ++k ) 
                 w2_inside[ k ] = 0; 
 
              //we use the transformation formula, therefor we have:
              for( int k = 0; k < dimension; ++k )
                for( int l = 0; l < dimension; ++l )
                  w2_inside[ k ] += (1.0 / epsilon_est_) * a_inside[ k ][ l ] * gradCor[ 0 ][ l ];


              //to save A_h(neighborEntityPoint) \nabla Corrector(Phi_H)(neighborEntityQuadPoint).
              RangeType w2_outside[dimension];
              for( int k = 0; k < dimension; ++k ) 
                 w2_outside[ k ] = 0; 
 
              //we use the transformation formula, therefor we have:
              for( int k = 0; k < dimension; ++k )
                for( int l = 0; l < dimension; ++l )
                  w2_outside[ k ] += (1.0 / epsilon_est_) * a_outside[ k ][ l ] * gradCorNeigh[ 0 ][ l ];


              //to save:
              //  A_h(neighborEntityPoint) [ \nabla Phi_H(globalPoint) 
              //                   + \nabla  Corrector(Phi_H)(neighborEntityQuadPoint) ].
              RangeType w_inside[dimension];
              for( int k = 0; k < dimension; ++k ) 
                 w_inside[ k ] = w1_inside[ k ] + w2_inside[ k ]; 

              //to save:
              //  A_h(neighborEntityPoint) [ \nabla Phi_H(globalPoint) 
              //                   + \nabla  Corrector(Phi_H)(neighborEntityQuadPoint) ].
              RangeType w_outside[dimension];
              for( int k = 0; k < dimension; ++k ) 
                 w_outside[ k ] = w1_outside[ k ] + w2_outside[ k ]; 

              //gradient jump:
              RangeType w_jump[dimension];
              for( int k = 0; k < dimension; ++k ) 
                 w_jump[ k ] = w_outside[ k ] - w_inside[ k ]; 

              //gradient jump * unitOuterNormal:
              RangeType w_unscaled = 0;
              for( int k = 0; k < dimension; ++k ) 
                 w_unscaled += w_jump[ k ] * unitOuterNormal[ k ]; 

              //gradient jump * scaledOuterNormal:
              RangeType w_scaled = 0;
              for( int k = 0; k < dimension; ++k ) 
                 w_scaled += w_jump[ k ] * scaledOuterNormal[ k ]; 
              // By scaling the second outer normal, we don't have to multiply the final result with the size of the face

              RangeType value;
              value = w_scaled * w_unscaled;
              value *= h_entity * h_entity * h_entity 
                       * innerFaceQuadrature.weight(faceQuadPoint);

              result += value;

             } // end loop face quadrature points

             } // end of if-condition:  nit.neighbor()
             else // (if nit is a boundary face:)
             { 
              //! hier verschenken wir nun etwas, solange wir den IntersectionIterator für nicht periodisches Gitter verwenden. Hier verlieren wir die 'Randsprünge', da wir 'zu früh' die L²-Norm bilden
              
              //to save gradient of Corrector(Phi_H) on entity:
              JacobianRangeType gradCor(0.0);
              //to save gradient of Corrector(Phi_H) on neighbor entity:
              JacobianRangeType gradCorNeigh(0.0);

              const PeriodicLocalFunctionType localFunction = corrector_u_H.localFunction(entity);  

              RangeType a[ dimension ][ dimension ];
              if (homogeneous_structure == true )
               {
                for( int k = 0; k < dimension; ++k )
                  for( int l = 0; l < dimension; ++l ) 
                    {
                      tensor_.evaluate( l, k,
                                        globalPoint,
                                        geometry.global( entityQuadrature.point( 0 ) ),
                                        a[ l ][ k ] );
                    }
               }
             else
                {

                 DomainType y_eps; // x_j + \eps_{est} * y
                 for( int k = 0; k < dimension; ++k )
                    y_eps[ k ] = globalPoint[ k ] + epsilon_est_ * geometry.global( entityQuadrature.point( 0 ) )[ k ];

                 for( int k = 0; k < dimension; ++k )
                   for( int l = 0; l < dimension; ++l ) 
                    {
                      tensor_.evaluate( l, k, y_eps, a[ l ][ k ] );
                    }
               }

             //to save A_h(EntityPoint) \nabla Phi_H(globalPoint):
             RangeType w1[dimension];
             for( int k = 0; k < dimension; ++k ) 
               w1[ k ] = 0;

             for( int k = 0; k < dimension; ++k )
               for( int l = 0; l < dimension; ++l )
                  w1[ k ] += a[ k ][ l ] * grad_u_H[ 0 ][ l ];


            FaceQuadratureType faceQuadrature(gridPart, nit, polOrd, FaceQuadratureType::INSIDE);

            const int quadratureNop = faceQuadrature.nop();
            for(int faceQuadPoint = 0; faceQuadPoint < quadratureNop; ++faceQuadPoint)
              {

               DomainType unitOuterNormal =
                 nit.unitOuterNormal(faceQuadrature.localPoint(faceQuadPoint));

              //outer normal Scaled with face size:
              DomainType scaledOuterNormal =
                nit.integrationOuterNormal(faceQuadrature.localPoint(faceQuadPoint));

              // get 'volume' of the face (this only works because we do not have curved faces):
              RangeType faceVolume(0.0);
              for ( int k = 0; k < dimension; ++k ) 
                faceVolume += scaledOuterNormal[k] * scaledOuterNormal[k];
              faceVolume = sqrt(faceVolume);

              // instead of using h_entity = entityVolume, we set:
              RangeType h_entity = entityVolume / faceVolume;
              h_entity = faceVolume; //! loeschen?
              // then it is independent of the dimension

              localFunction.jacobian( faceQuadrature , faceQuadPoint , gradCor );

              //to save A_h(EntityPoint) \nabla Corrector(Phi_H)(EntityQuadPoint).
              RangeType w2[dimension];
              for( int k = 0; k < dimension; ++k ) 
                 w2[ k ] = 0; 
 
              //we use the transformation formula, therefor we have:
              for( int k = 0; k < dimension; ++k )
                for( int l = 0; l < dimension; ++l )
                  w2[ k ] += (1.0 / epsilon_est_) * a[ k ][ l ] * gradCor[ 0 ][ l ];

              //to save:
              //  A_h(neighborEntityPoint) [ \nabla Phi_H(globalPoint) 
              //                   + \nabla  Corrector(Phi_H)(neighborEntityQuadPoint) ].
              RangeType w[dimension];
              for( int k = 0; k < dimension; ++k ) 
                 w[ k ] = w1[ k ] + w2[ k ]; 


              // * unitOuterNormal:
              RangeType w_unscaled = 0;
              for( int k = 0; k < dimension; ++k ) 
                 w_unscaled += w[ k ] * unitOuterNormal[ k ]; 

              //gradient jump * scaledOuterNormal:
              RangeType w_scaled = 0;
              for( int k = 0; k < dimension; ++k ) 
                 w_scaled += w[ k ] * scaledOuterNormal[ k ]; 
              // By scaling the second outer normal, we don't have to multiply the final result with the size of the face

              RangeType value;
              value = w_scaled * w_unscaled;
              value *= h_entity * h_entity * h_entity 
                       * faceQuadrature.weight(faceQuadPoint);

              result += value;

             }

            } // 'end of else'

            } // end of intersection iterator loop

      } // end of entity iterator loop

      // each face was visited twice
      result = result / 2.0; 

      return result;

    } // end of method 


    // return local approximation error (already locally squared)
    RangeType errorMu( const EntityType &entity,
                       const PeriodicDiscreteFunctionType &corrector_u_H,
                       const JacobianRangeType &gradient_u_H_in_xj )
     {


      bool homogeneous_structure = tensor_.homogeneous();

      // by global quadrature we mean: a quadrature on the macro element T_j,
      // by local quadrature we mean: a quadrature on the unit cube Y

      EntityQuadratureType globalQuadrature(entity, 0 ); // 0 = polynomial order
      //for further considerations with matrices of the type A(x,y) (whereas the dependency on the first variable x is a 'real' dependency!) it may be better to use polynomial order = 2 instead of polynomial order = 0!
      // for matrices of type A(x,y) = A(y), it is of course unnecessary to use polynomial order = 2

      const EntityGeometryType& globalGeometry = entity.geometry();

      const DomainType &xj = globalGeometry.global(globalQuadrature.point(0));

      RangeType errorMu(0.);

      const int quadratureNop = globalQuadrature.nop();
      for(int quadPointIn_T_j = 0; quadPointIn_T_j < quadratureNop; ++quadPointIn_T_j)
      // start loop over the quadrature points of the macro mesh element T_j
       {

        // to save the approximation of
        //     \int_Y ((A - A_h)(\nabla_x u_H + \nabla_y K_h(u_H)))^2
        // in corresponding quadrature points:
        RangeType valueOfLocalQuadrature(0.);

        PeriodicIteratorType endit = periodicDiscreteFunctionSpace_.end();
        for(PeriodicIteratorType it = periodicDiscreteFunctionSpace_.begin(); it != endit ; ++it)
        // start loop over the elements of the unit cell mesh
         {

           // entity
           const PeriodicEntityType& localEntity = *it;

           // create quadrature for given geometry type 
           PeriodicEntityQuadratureType localQuadrature(localEntity , 2 ); // 2 = polynomial order
           // create a trivial quadrature
           // This is important to compare A and A_h properly. A_h is defined by the values of A in the quadrature points of the trivial quadrature: A_h(x,y)_{|T_j /times S_k^*} := A(x_j, y_k). Now if we used the same quadrature for A, then the error would be zero! (A_h(x_j, y_k) - A(x_j, y_k) = 0) This is wrong since the error should have the same order as theory suggests. (generaly first order)
           PeriodicEntityQuadratureType trivialQuadrature(localEntity , 0 ); // 0 = polynomial order

           // get geoemetry of local entity
           const PeriodicEntityGeometryType& localGeometry = localEntity.geometry();

           PeriodicLocalFunctionType cor_u_H = corrector_u_H.localFunction(localEntity); 

           // save A_h(x_j, y_k)
           RangeType a_h[ dimension ][ dimension]; 
           if (homogeneous_structure == true )
            {
             // evaluate function in global quadrature Point:
             for( int k = 0; k < dimension ; ++k )
               for( int l = 0; l < dimension ; ++l ) 
                 tensor_.evaluate( k, l, xj,
                                   localGeometry.global( trivialQuadrature.point( 0 ) ),
                                   //0, since we use the trivial quadrature
                                   a_h[ k ][ l ] );
            }
           else
            {
              DomainType y_eps; // x_j + \eps_{est} * y
              for( int k = 0; k < dimension; ++k )
                 y_eps[ k ] = xj[ k ] + epsilon_est_ * localGeometry.global( trivialQuadrature.point( 0 ) )[ k ];

              // evaluate function in global quadrature Point:
              for( int k = 0; k < dimension ; ++k )
               for( int l = 0; l < dimension ; ++l ) 
                 tensor_.evaluate( k, l, y_eps, a_h[ k ][ l ] );
            }

           // we need to average over the local distances, for details see end of localQuadPoint-loop
           double average = 0;

           // integrate 
           const int localQuadratureNop = localQuadrature.nop();
           for(int localQuadPoint = 0; localQuadPoint < localQuadratureNop; ++localQuadPoint)
            {

              // save A(x,y) in corresponding, non trivial, quadrature points
              RangeType a[ dimension ][ dimension ];
              if (homogeneous_structure == true )
               {
                for( int k = 0; k < dimension ; ++k )
                  for( int l = 0; l < dimension ; ++l ) 
                    tensor_.evaluate( k, l,
                                      globalGeometry.global( globalQuadrature.point( quadPointIn_T_j ) ),
                                      localGeometry.global( localQuadrature.point( localQuadPoint ) ),
                                      a[ k ][ l ] );
               }
              else
               {
                 DomainType y_eps_new;
                 for( int k = 0; k < dimension; ++k )
                 y_eps_new[ k ] = globalGeometry.global( globalQuadrature.point( quadPointIn_T_j ) )[ k ] + epsilon_est_ * localGeometry.global( localQuadrature.point( localQuadPoint ) )[ k ];

                 for( int k = 0; k < dimension ; ++k )
                  for( int l = 0; l < dimension ; ++l ) 
                    tensor_.evaluate( k, l, y_eps_new, a[ k ][ l ] );
               }

              PeriodicJacobianRangeType gradient_Cor_u_H;
              cor_u_H.jacobian( localQuadrature , localQuadPoint , gradient_Cor_u_H );

              const double locQuadWeight = localQuadrature.weight(localQuadPoint) *
                  localGeometry.integrationElement(localQuadrature.point(localQuadPoint));

              RangeType approxerror[ dimension ];
              for( int k = 0; k < dimension; ++k )
                 approxerror[ k ] = 0;


              for( int k = 0; k < dimension; ++k )
               for( int l = 0; l < dimension; ++l )
                {
                 approxerror[ k ] += ( a[ k ][ l ] - a_h[ k ][ l ] )
                   * ( gradient_u_H_in_xj[ 0 ][ l ] + (1.0 / epsilon_est_) * gradient_Cor_u_H[ 0 ][ l ] ); 
                }


              RangeType val(0);
              for( int k = 0; k < dimension; ++k )
               val += approxerror[ k ] * approxerror[ k ];

              average += val;

             valueOfLocalQuadrature += locQuadWeight;

            } // end of localQuadPoint-loop

            // averaging prevents that the maximum of all distances (A - A_h) is multiplied with the maximum of the weights, which could enlarge the error and decrease the rate of convergence fundamentaly!
            average = average / localQuadratureNop;
            valueOfLocalQuadrature *= average;

          } // end loop over the elements of the unit cell mesh

          const double weight = globalQuadrature.weight(quadPointIn_T_j) *
             globalGeometry.integrationElement(globalQuadrature.point(quadPointIn_T_j));

          errorMu += weight * valueOfLocalQuadrature;

        } // end loop over the quadrature points of the macro mesh element T_j

     return errorMu;

     } // end of method 


    RangeType errorZeta( const EntityType &entity,
                               DiscreteFunctionType &u_H,
                         const PeriodicDiscreteFunctionType &correctorU_H_onEntity,
                         const JacobianRangeType &gradientU_H_onEntity )
     {
      RangeType errorZeta(0.);

      //const EntityType& entity = *entityPointer;

      bool homogeneous_structure = tensor_.homogeneous();

      double epsilon_est;
      tensor_.getEpsilonEstimated(epsilon_est);

      EntityQuadratureType entityQuadrature(entity,0); // 0 = polynomial order
      // the global quadrature (quadrature on the macro element T_j)

      const EntityGeometryType& globalEntityGeometry = entity.geometry();

      const DomainType &xj = globalEntityGeometry.global(entityQuadrature.point(0));

      const RangeType entityVolume = entityQuadrature.weight(0) * 
              globalEntityGeometry.integrationElement(entityQuadrature.point(0));

      double cof[10000];
      for( int r = 0; r < 10000; ++r )
         {cof[r] = 0;} 

      int counting = 0;
      // we save the value of the hmfem solution in its dofs.
      const DofIteratorType end = u_H.dend();
      for( DofIteratorType dit = u_H.dbegin(); dit != end; ++dit )
        { cof[ counting ] = *dit;
          counting += 1;
        }

      bool reader_is_open = false;

      DiscreteFunctionReader dfr( (cell_solution_location_).c_str() );

      reader_is_open = dfr.open();

      // to add and multiply discfunctions
      DiscFuncAdapter<PeriodicDiscreteFunctionType> adapter;



//!alternative (testen)! Vergleich mit dem uebergebenen correctorU_H_onEntity!
#if 0
         const BaseFunctionSetType baseSet
            = discreteFunctionSpace_.baseFunctionSet( entity );
         const int numBaseFunctions = baseSet.numBaseFunctions();

         // we want to determine correctorU_H:
         PeriodicDiscreteFunctionType correctorU_H_onEntity( "correctorU_H" , periodicDiscreteFunctionSpace_);
         correctorU_H_onEntity.clear();

         // an auxilliary function (reconstruction of a certain base function), that is just used for saving certain values
         PeriodicDiscreteFunctionType correctorPhi_i_onEntity( "corrector Phi_i" , periodicDiscreteFunctionSpace_ );

         // for-loop to determine correctorU_H (on actual entity)
         for( int i = 0; i < numBaseFunctions; ++i )
          {
            int number_of_cell_prob = get_cell_problem_id( entityPointer, i );

            if (reader_is_open)
              {dfr.read( number_of_cell_prob, correctorPhi_i_onEntity );}

            int globalNumberOfBaseFunction = discreteFunctionSpace_.mapToGlobal( entity , i );

//! z.B. coefficient[ globalNumberOfBaseFunction ] -> - coefficient[ globalNumberOfBaseFunction ]
//! dann beide Varianten addieren und es muss Null rauskommen!

            // K_h(PHI_i) = u_i * K_h(PHI_i)
            adapter.multScalar( correctorPhi_i_onEntity , cof[ globalNumberOfBaseFunction ] );
            // K_h(u_H) = K_h(u_H) + K_h(PHI_i)
            adapter.add( correctorU_H_onEntity , correctorPhi_i_onEntity );

          } // end for loop
#endif



      const GridPartType &gridPart = discreteFunctionSpace_.gridPart();


      PeriodicDiscreteFunctionType
                correctorPhi_i_onNeighEntity( "correctorPhi_i_onNeighEntity",
                                              periodicDiscreteFunctionSpace_ );

      PeriodicDiscreteFunctionType
             correctorU_H_onNeighEntity( "correctorU_H_onNeighEntity" ,
                                         periodicDiscreteFunctionSpace_ );

      enum { maxnumOfBaseFct = 100 };

      JacobianRangeType gradientU_H_onNeighEntity(0.);

      IntersectionIteratorType endnit = gridPart.iend(entity);
      for(IntersectionIteratorType nit = gridPart.ibegin(entity); nit != endnit ; ++nit)
         {
          if ( nit.neighbor() ) //if there is a neighbor entity
           {
             correctorU_H_onNeighEntity.clear();
             for( int k = 0; k < dimension; ++k )
                gradientU_H_onNeighEntity[k] = 0;

             const EntityPointerType neighborEntityPointer = nit.outside();
             const EntityType& neighborEntity = *neighborEntityPointer;
             // get quadrature of neighbor entity
             EntityQuadratureType neighborEntityQuadrature(neighborEntity,0); //PolOrd = 0 
             // only one qudraturePoint!
             const EntityGeometryType& globalNeighEntityGeometry = neighborEntity.geometry();
 
             LocalFunctionType u_H_neigh_entity = u_H.localFunction(neighborEntity); 

             const RangeType neighEntityVolume = neighborEntityQuadrature.weight(0) * 
                 globalNeighEntityGeometry.integrationElement(neighborEntityQuadrature.point(0));

             const BaseFunctionSetType neighEntityBaseSet
                 = discreteFunctionSpace_.baseFunctionSet( neighborEntity );
             const int numBaseFunctionsOnEntityNeigh = neighEntityBaseSet.numBaseFunctions();

//! Start computing the (unit) outer normal belonging to the current face (of the global mesh!)
             FaceQuadratureType faceQuadrature(gridPart, nit, 0, FaceQuadratureType::INSIDE);

             // unit outer normal to the global entity face
             DomainType unitOuterNormal =
                nit.unitOuterNormal(faceQuadrature.localPoint( 0 ));

              //outer normal Scaled with face size:
             DomainType scaledOuterNormal =
                nit.integrationOuterNormal(faceQuadrature.localPoint( 0 ));

             RangeType face_size = 0.0;
             for( int l = 0; l < dimension; ++l ) 
              {
               if (fabs(scaledOuterNormal[l]) >= face_size )
                 { face_size = fabs(scaledOuterNormal[l]); }
              }

//! anpassen: localGeometry.integrationElement(localQuadrature.point( 0 ));

//! End computing outer normal

             u_H_neigh_entity.jacobian( neighborEntityQuadrature , 0 , gradientU_H_onNeighEntity );

             for( int i = 0; i < numBaseFunctionsOnEntityNeigh; ++i )
              {

                int globalNumberOfBaseFunction = discreteFunctionSpace_.mapToGlobal( neighborEntity , i );

                int number_of_cell_prob = get_cell_problem_id( neighborEntityPointer, i );

                correctorPhi_i_onNeighEntity.clear();
                if (reader_is_open)
                  {dfr.read( number_of_cell_prob, correctorPhi_i_onNeighEntity );}

                // K_h(PHI_i) = u_i * K_h(PHI_i)
                adapter.multScalar( correctorPhi_i_onNeighEntity , cof[ globalNumberOfBaseFunction ] );
                // K_h(u_H) = K_h(u_H) + K_h(PHI_i)
                adapter.add( correctorU_H_onNeighEntity , correctorPhi_i_onNeighEntity );

              } 

             // Now gradientU_H_onNeighEntity and correctorU_H_onNeighEntity are determined, the corresponding base functions on neighborEntity are no more important now.


             RangeType averaged_jump = 0.0;

             // Start
             PeriodicIteratorType endit = periodicDiscreteFunctionSpace_.end();
             for(PeriodicIteratorType it = periodicDiscreteFunctionSpace_.begin(); it != endit ; ++it)
               // start loop over the elements of the unit cell mesh
               {
                 const PeriodicEntityType& localEntity = *it;
                 // by local entity we mean an entity element of the local mesh on the unit cell Y

                 // create quadrature for given geometry type 
                 PeriodicEntityQuadratureType localQuadrature(localEntity , 0 ); //polOrd = 0

                 // get geoemetry of local entity
                 const PeriodicEntityGeometryType& localGeometry = localEntity.geometry();


                 PeriodicLocalFunctionType cor_u_H_entity
                      = correctorU_H_onEntity.localFunction(localEntity);

                 PeriodicJacobianRangeType grad_cor_u_H_entity;
                 cor_u_H_entity.jacobian
                     ( localQuadrature , 0 /*localQuadPoint*/ , grad_cor_u_H_entity );


                 PeriodicLocalFunctionType cor_u_H_neighEntity
                      = correctorU_H_onNeighEntity.localFunction(localEntity);

                 PeriodicJacobianRangeType grad_cor_u_H_neighEntity;
                 cor_u_H_neighEntity.jacobian
                     ( localQuadrature , 0 /*localQuadPoint*/ , grad_cor_u_H_neighEntity );

 
//! Start computing A_h( x_entity, y_k ) and A_h( x_neighEntity, y_k )

// y_k is the current (local) quadrature point
                 RangeType a_entity[ dimension ][ dimension ];
                 if (homogeneous_structure == true )
                  {
                    for( int k = 0; k < dimension; ++k )
                     for( int l = 0; l < dimension; ++l ) 
                      {
                       tensor_.evaluate( l, k, xj, //global quad point
                                         localGeometry.global( localQuadrature.point( 0 ) ), //local quad point
                                         a_entity[ l ][ k ] );
                      }
                  }
                 else
                  {

                    DomainType y_eps_entity; // x_j + \eps_{est} * y
                    for( int k = 0; k < dimension; ++k )
                      y_eps_entity[ k ] = xj[ k ] + epsilon_est_ * localGeometry.global( localQuadrature.point( 0 ) )[ k ];

                    for( int k = 0; k < dimension; ++k )
                     for( int l = 0; l < dimension; ++l ) 
                      {
                       tensor_.evaluate( l, k, y_eps_entity,
                                         a_entity[ l ][ k ] );
                      }
                  }

                 RangeType a_neighEntity[ dimension ][ dimension ];
                 if (homogeneous_structure == true )
                  {
                    for( int k = 0; k < dimension; ++k )
                     for( int l = 0; l < dimension; ++l ) 
                      {
                       tensor_.evaluate( l, k,
                                         globalNeighEntityGeometry.global( neighborEntityQuadrature.point(0) ), //global quad point
                                          localGeometry.global( localQuadrature.point( 0 ) ), //local quad point
                                          a_neighEntity[ l ][ k ] );
                      }
                  }
                 else
                  {
                    DomainType y_eps_neigh_entity; // x_j + \eps_{est} * y
                    for( int k = 0; k < dimension; ++k )
                      y_eps_neigh_entity[ k ] = globalNeighEntityGeometry.global( neighborEntityQuadrature.point(0) )[ k ] + epsilon_est_ * localGeometry.global( localQuadrature.point( 0 ) )[ k ];


                    for( int k = 0; k < dimension; ++k )
                     for( int l = 0; l < dimension; ++l ) 
                      {
                       tensor_.evaluate( l, k, y_eps_neigh_entity,
                                          a_neighEntity[ l ][ k ] );
                      }
                  }
//! End computing A_h( x_entity, y_k ) and A_h( x_neighEntity, y_k )

//! Start computing the gradient jump:
//!  [ A_h( x_entity, y_k ) ( \nabla u_H( x_entity ) + \nabla K_h(u_H) ( x_entity, y_k ) )
//!  - A_h( x_neighEntity, y_k ) ( \nabla u_H( x_neighEntity ) + \nabla K_h(u_H) ( x_neighEntity, y_k ) ) ]
//!        * n_entity

                 //to save  A_h( x_entity, y_k ) \nabla u_H( x_entity ):
                 RangeType A_h_grad_u_H_entity[dimension];
                 for( int k = 0; k < dimension; ++k ) 
                  A_h_grad_u_H_entity[ k ] = 0; 

                 for( int k = 0; k < dimension; ++k )
                  for( int l = 0; l < dimension; ++l )
                   A_h_grad_u_H_entity[ k ] += a_entity[ k ][ l ] * gradientU_H_onEntity[ 0 ][ l ];



                 //to save  A_h( x_neighEntity, y_k ) \nabla u_H( x_neighEntity ):
                 RangeType A_h_grad_u_H_neighEntity[dimension];
                 for( int k = 0; k < dimension; ++k ) 
                  A_h_grad_u_H_neighEntity[ k ] = 0; 

                 for( int k = 0; k < dimension; ++k )
                  for( int l = 0; l < dimension; ++l )
                   A_h_grad_u_H_neighEntity[ k ] += a_neighEntity[ k ][ l ] * gradientU_H_onNeighEntity[ 0 ][ l ];


                 //to save A_h( x_entity, y_k ) \nabla K_h(u_H)( x_entity, y_k )
                 RangeType A_h_grad_COR_u_H_entity[dimension];
                 for( int k = 0; k < dimension; ++k ) 
                   A_h_grad_COR_u_H_entity[ k ] = 0; 

                 //since we did not really compute K_h(u_H) but some K_h(u_H)(x_j + \epsilon y)
                 //we need to use the transformation formula. Therefor we have:
                 for( int k = 0; k < dimension; ++k )
                   for( int l = 0; l < dimension; ++l )
                     A_h_grad_COR_u_H_entity[ k ] +=
                      (1.0 / epsilon_est) * a_entity[ k ][ l ] * grad_cor_u_H_entity[ 0 ][ l ];


                 //to save A_h( x_neighEntity, y_k ) \nabla K_h(u_H)( x_neighEntity, y_k )
                 RangeType A_h_grad_COR_u_H_neighEntity[dimension];
                 for( int k = 0; k < dimension; ++k ) 
                   A_h_grad_COR_u_H_neighEntity[ k ] = 0; 

                 for( int k = 0; k < dimension; ++k )
                   for( int l = 0; l < dimension; ++l )
                     A_h_grad_COR_u_H_neighEntity[ k ] +=
                      (1.0 / epsilon_est) * a_neighEntity[ k ][ l ] * grad_cor_u_H_neighEntity[ 0 ][ l ];


                 //to save:
                 //  A_h( x_entity, y_k ) ( \nabla u_H( x_entity ) + \nabla K_h(u_H) ( x_entity, y_k )
                 RangeType w_entity[dimension];
                 for( int k = 0; k < dimension; ++k ) 
                   w_entity[ k ] = A_h_grad_u_H_entity[ k ] + A_h_grad_COR_u_H_entity[ k ]; 
                 // REC stands for Reconstruction 

                 //to save:
                 //  A_h( x_neighEntity, y_k ) [ \nabla u_H( x_neighEntity ) 
                 //            + \nabla K_h(u_H) ( x_neighEntity, y_k ) ]
                 RangeType w_neighEntity[dimension];
                 for( int k = 0; k < dimension; ++k ) 
                   w_neighEntity[ k ] = A_h_grad_u_H_neighEntity[ k ] + A_h_grad_COR_u_H_neighEntity[ k ]; 

                 //to save: jump * 'size of the global face'
                 RangeType jump(0.);
                 for( int k = 0; k < dimension; ++k )
                   jump += (w_entity[ k ] - w_neighEntity[ k ])
                            * unitOuterNormal[ k ];

//! End computing gradient jump

                 const double localweight = localQuadrature.weight( 0 ) *
                   localGeometry.integrationElement(localQuadrature.point( 0 ));

                 averaged_jump += jump * localweight;

               }

             errorZeta += averaged_jump * averaged_jump * faceQuadrature.weight( 0 ) * face_size;

           } // end of if-loop (if there is a neighbor entity...)

         } //end global intersection Iterator

     // every gradient-jump will occur twice:
     errorZeta = errorZeta / 2.0;

     return errorZeta;

     } // end of method

}; // end of class ErrorEstimater

} // end namespace 

#endif
