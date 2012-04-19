#ifndef DUNE_MSFEM_ERRORESTIMATOR_HH
#define DUNE_MSFEM_ERRORESTIMATOR_HH

//#include <dune/common/classname.hh>

// where the quadratures are defined 
#include <dune/fem/quadrature/cachingquadrature.hh>

#include <dune/multiscale/tools/errorestimation/MsFEM/conservative_flux_solver.hh>


namespace Dune 
{


  template < class DiscreteFunctionImp,
             class DiffusionImp,
             class SourceImp, 
             class MacroMicroGridSpecifierImp,
             class SubGridListImp >
  class MsFEMErrorEstimator
  {

    typedef DiffusionImp DiffusionOperatorType;
    typedef SourceImp SourceType;
    typedef MacroMicroGridSpecifierImp MacroMicroGridSpecifierType;
    typedef SubGridListImp SubGridListType;
    
    //! Necessary typedefs for the DiscreteFunctionImp:

    typedef DiscreteFunctionImp DiscreteFunctionType;
    typedef typename DiscreteFunctionType :: FunctionSpaceType FunctionSpaceType;

    typedef typename DiscreteFunctionType :: LocalFunctionType LocalFunctionType;
    typedef typename DiscreteFunctionType :: DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

    typedef typename DiscreteFunctionSpaceType :: RangeFieldType RangeFieldType;
    typedef typename DiscreteFunctionSpaceType :: DomainType DomainType;
    typedef typename DiscreteFunctionSpaceType :: RangeType RangeType;
    typedef typename DiscreteFunctionSpaceType :: JacobianRangeType
      JacobianRangeType;

    typedef typename DiscreteFunctionType :: DofIteratorType DofIteratorType;
    typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType;
    typedef typename DiscreteFunctionSpaceType :: GridType GridType;
    typedef typename DiscreteFunctionSpaceType :: IteratorType IteratorType;
    typedef typename DiscreteFunctionSpaceType :: LagrangePointSetType LagrangePointSetType;
    
    typedef typename GridType :: Traits :: LeafIndexSet LeafIndexSetType;
    
    typedef typename GridPartType :: IntersectionIteratorType IntersectionIteratorType;
    typedef typename IntersectionIteratorType::Intersection Intersection;
    typedef typename GridType :: template Codim<0> :: Entity EntityType; 
    typedef typename GridType :: template Codim<0> :: EntityPointer            EntityPointerType;
    typedef typename GridType :: template Codim<1> :: Entity FaceType; 
    typedef typename GridType :: template Codim<1> :: EntityPointer            FacePointerType;    
    typedef typename GridType :: template Codim<0> :: Geometry                 EntityGeometryType; 
    typedef typename GridType :: template Codim<1> :: Geometry                 FaceGeometryType;
    typedef typename DiscreteFunctionSpaceType     :: BaseFunctionSetType      BaseFunctionSetType;
    
    typedef CachingQuadrature < GridPartType , 0 > EntityQuadratureType;
    typedef CachingQuadrature < GridPartType , 1 > FaceQuadratureType;


    // --------------------------- subgrid typedefs ------------------------------------

    typedef SubGrid< GridType::dimension , GridType > SubGridType; 
    typedef LeafGridPart< SubGridType > SubGridPart;

    typedef typename SubGridType :: Traits :: LeafIndexSet SubGridLeafIndexSet;

    typedef LagrangeDiscreteFunctionSpace < FunctionSpaceType, SubGridPart, 1 > //1=POLORDER
      SubGridDiscreteFunctionSpaceType;

   typedef AdaptiveDiscreteFunction < SubGridDiscreteFunctionSpaceType > SubGridDiscreteFunctionType;

   typedef typename SubGridDiscreteFunctionSpaceType :: IteratorType SubGridIteratorType;

   typedef typename SubGridIteratorType :: Entity SubGridEntityType;

   typedef typename SubGridEntityType :: EntityPointer SubGridEntityPointerType; 

   typedef typename SubGridDiscreteFunctionType :: LocalFunctionType SubGridLocalFunctionType;

   typedef typename SubGridDiscreteFunctionSpaceType :: LagrangePointSetType
    SubGridLagrangePointSetType;

   enum { faceCodim = 1 };
   typedef typename SubGridLagrangePointSetType :: template Codim< faceCodim > 
                                                :: SubEntityIteratorType
   SubGridFaceDofIteratorType;
    
    
   typedef typename LagrangePointSetType :: template Codim< faceCodim > 
                                         :: SubEntityIteratorType
   FaceDofIteratorType;
    

   typedef typename SubGridType :: template Codim<0> :: Geometry                 SubGridEntityGeometryType; 

   //!-----------------------------------------------------------------------------------------


    enum { dimension = GridType :: dimension};
    enum { spacePolOrd = DiscreteFunctionSpaceType :: polynomialOrder }; 
    enum { maxnumOfBaseFct = 100 };

#if 0


  private:

    const PeriodicDiscreteFunctionSpaceType &periodicDiscreteFunctionSpace_;
    const DiscreteFunctionSpaceType &discreteFunctionSpace_;
    const DiscreteFunctionSpaceType &auxiliaryDiscreteFunctionSpace_;
    // an auxiliaryDiscreteFunctionSpace to get an Intersection Iterator for the periodicDiscreteFunctionSpace
    // (for the periodic grid partition there is no usable intersection iterator implemented, therefore we use the intersection iterator for the corresponding non-periodic grid partition (this not efficient and increases the estimated error, but it works)
#endif

    const DiscreteFunctionSpaceType& fineDiscreteFunctionSpace_;
    MacroMicroGridSpecifierType& specifier_;
    SubGridListType& subgrid_list_;
    const DiffusionOperatorType &diffusion_;
    const SourceType& f_;
    std :: string& path_;


  public:

    MsFEMErrorEstimator ( const DiscreteFunctionSpaceType &fineDiscreteFunctionSpace,
                          MacroMicroGridSpecifierType& specifier,
                          SubGridListType& subgrid_list,
                          const DiffusionOperatorType& diffusion,
                          const SourceType& f,
                          std :: string& path )
    : fineDiscreteFunctionSpace_( fineDiscreteFunctionSpace ),
      specifier_( specifier ),
      subgrid_list_( subgrid_list ),
      diffusion_( diffusion ),
      f_( f ),
      path_( path )
    {
    }


    // create a hostgrid function from a subgridfunction
    void subgrid_to_hostrid_function( const SubGridDiscreteFunctionType &sub_func,
                                            DiscreteFunctionType &host_func )
    {

       host_func.clear();

       const SubGridDiscreteFunctionSpaceType &subDiscreteFunctionSpace = sub_func.space();
       const SubGridType &subGrid = subDiscreteFunctionSpace.grid();       
       
       SubGridIteratorType sub_endit = subDiscreteFunctionSpace.end();
       for( SubGridIteratorType sub_it = subDiscreteFunctionSpace.begin(); sub_it != sub_endit; ++sub_it )
          {

             const SubGridEntityType &sub_entity = *sub_it;

             EntityPointerType host_entity_pointer = subGrid.template getHostEntity<0>( *sub_it );
             const EntityType& host_entity = *host_entity_pointer;

             SubGridLocalFunctionType sub_loc_value = sub_func.localFunction( sub_entity );
             LocalFunctionType host_loc_value = host_func.localFunction( host_entity );

             const unsigned int numBaseFunctions = sub_loc_value.baseFunctionSet().numBaseFunctions();
             for( unsigned int i = 0; i < numBaseFunctions; ++i )
               {
                 host_loc_value[ i ] = sub_loc_value[ i ];
               }

          }
    }





    // create twoa hostgrid functions from two subgridfunctions
    void subgrid_to_hostrid_function( const SubGridDiscreteFunctionType &sub_func_1,
                                      const SubGridDiscreteFunctionType &sub_func_2,
                                            DiscreteFunctionType &host_func_1,
                                            DiscreteFunctionType &host_func_2 )
    {

       host_func_1.clear();
       host_func_2.clear();

       const SubGridDiscreteFunctionSpaceType &subDiscreteFunctionSpace = sub_func_1.space();
       const SubGridType &subGrid = subDiscreteFunctionSpace.grid();       
       
       SubGridIteratorType sub_endit = subDiscreteFunctionSpace.end();
       for( SubGridIteratorType sub_it = subDiscreteFunctionSpace.begin(); sub_it != sub_endit; ++sub_it )
          {

             const SubGridEntityType &sub_entity = *sub_it;

             EntityPointerType host_entity_pointer = subGrid.template getHostEntity<0>( *sub_it );
             const EntityType& host_entity = *host_entity_pointer;

             SubGridLocalFunctionType sub_loc_value_1 = sub_func_1.localFunction( sub_entity );
             SubGridLocalFunctionType sub_loc_value_2 = sub_func_2.localFunction( sub_entity );
             LocalFunctionType host_loc_value_1 = host_func_1.localFunction( host_entity );
             LocalFunctionType host_loc_value_2 = host_func_2.localFunction( host_entity );

             const unsigned int numBaseFunctions = sub_loc_value_1.baseFunctionSet().numBaseFunctions();
             for( unsigned int i = 0; i < numBaseFunctions; ++i )
               {
                 host_loc_value_1[ i ] = sub_loc_value_1[ i ];
                 host_loc_value_2[ i ] = sub_loc_value_2[ i ];
               }

          }
    }




    //! method to get the local mesh size H of a coarse grid entity 'T'
    // works only for our 2D examples!!!! 
    RangeType get_coarse_grid_H( const EntityType &entity)
    {

      // entity_H means H (the diameter of the entity)
      RangeType entity_H = 0.0;

      const GridPartType &coarseGridPart = specifier_.coarseSpace().gridPart();

      // compute the size of the faces of the entities and selected the largest.
      IntersectionIteratorType endnit = coarseGridPart.iend(entity);
      for( IntersectionIteratorType nit = coarseGridPart.ibegin(entity); nit != endnit ; ++nit)
        {
          FaceQuadratureType innerFaceQuadrature( coarseGridPart, *nit, 0 , FaceQuadratureType::INSIDE);

          DomainType scaledOuterNormal =
               nit->integrationOuterNormal(innerFaceQuadrature.localPoint(0));

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



    // for a coarse grid entity T:
    // return:  H_T ||f||_{L^2(T)}
    RangeType indicator_f( const EntityType &entity )
    {

        // create quadrature for given geometry type 
        CachingQuadrature <GridPartType , 0 > entityQuadrature(entity, 2*spacePolOrd+2 ); 

        // get geoemetry of entity
        const EntityGeometryType& geometry = entity.geometry();

        RangeType H_T = get_coarse_grid_H(entity);

        RangeType y(0);
        RangeType local_indicator(0);

        const int quadratureNop = entityQuadrature.nop();
        for(int quadraturePoint = 0; quadraturePoint < quadratureNop; ++quadraturePoint)
         {
          const double weight = entityQuadrature.weight(quadraturePoint) * 
              geometry.integrationElement(entityQuadrature.point(quadraturePoint));

          f_.evaluate( geometry.global( entityQuadrature.point(quadraturePoint) ),y );
          y = y * y;

          local_indicator += weight * y;
         }

        local_indicator *= pow(H_T, 2.0);

        return sqrt(local_indicator);
    }

    // is a given point on a given face?
    bool point_on_face( const Intersection& face, const DomainType& point )
     {

       DomainType corner_0 = face.geometry().corner(0);
       DomainType corner_1 = face.geometry().corner(1);

       if ( corner_0[0] == corner_1[0] )
        {
          if ( point[0] != corner_0[0] )
            { return false; }
          else
            {
              RangeType lambda = ( point[1] - corner_1[1] ) / ( corner_0[1] - corner_1[1] );
              if ( (lambda >= 0.0) && (lambda <= 1.0) )
               { return true; }
              else
               { return false; }
            }
        }
       else
        {
          RangeType lambda = ( point[0] - corner_1[0] ) / ( corner_0[0] - corner_1[0] );

          if ( (lambda >= 0.0) && (lambda <= 1.0) )
           {
             RangeType convex_comb = (lambda * corner_0[1]) + ((1.0 - lambda) * corner_1[1]);
             if ( convex_comb == point[1] )
              { return true; }
             else
              { return false; }
           }
          else
           { return false; }

        }
     }

    // is a given face part of another given coarse face
    bool is_subface( const Intersection& fine_face, const Intersection& coarse_face )
     {

       DomainType corner_0 = fine_face.geometry().corner(0);
       DomainType corner_1 = fine_face.geometry().corner(1);

       if ( point_on_face( coarse_face, corner_0 ) && point_on_face( coarse_face, corner_1 ) )
        { return true; }
       else
        { return false; }

     }

#if 1
    // jump in conservative flux
    void getFluxes( const EntityType &coarse_entity,
                    const DiscreteFunctionType& msfem_coarse_part,
                          RangeType& jump_conservative_flux,
                          RangeType& jump_coarse_flux )
     {
       
       // jump for each face
       RangeType jump[3];

       jump[0] = 0.0;
       jump[1] = 0.0;
       jump[2] = 0.0;
       
       // coarse grid jump for each face
       RangeType coarse_jump[3];
       
       coarse_jump[0] = 0.0;
       coarse_jump[1] = 0.0;
       coarse_jump[2] = 0.0;
       
       const DiscreteFunctionSpaceType& coarseDiscreteFunctionSpace = specifier_.coarseSpace();
       const LeafIndexSetType& coarseGridLeafIndexSet = coarseDiscreteFunctionSpace.gridPart().grid().leafIndexSet();

       int index_coarse_entity = coarseGridLeafIndexSet.index( coarse_entity );

       //! ---- get sub grid that for coarse entity

       // the sub grid U(T) that belongs to the coarse_grid_entity T
       SubGridType& sub_grid_U_T = subgrid_list_.getSubGrid( index_coarse_entity );
       SubGridPart subGridPart( sub_grid_U_T );

       SubGridDiscreteFunctionSpaceType localDiscreteFunctionSpace( subGridPart );

       SubGridDiscreteFunctionType conservative_flux_coarse_ent_e0( "Conservative Flux on coarse entity for e_0", localDiscreteFunctionSpace );
       conservative_flux_coarse_ent_e0.clear();

       SubGridDiscreteFunctionType conservative_flux_coarse_ent_e1( "Conservative Flux on coarse entity for e_1", localDiscreteFunctionSpace );
       conservative_flux_coarse_ent_e1.clear();

       // --------- load local solutions -------

       char location_lps_0[50];
       sprintf( location_lps_0, "_conservativeFlux_e_%d_sg_%d", 0, index_coarse_entity );
       std::string location_lps_s_0( location_lps_0 );

       std :: string cf_solution_location_0;

       // the file/place, where we saved the solutions conservative flux problems problems
       cf_solution_location_0 = path_ + "/cf_problems/" + location_lps_s_0;

       bool reader_is_open = false;
       // reader for data file:
       DiscreteFunctionReader discrete_function_reader_0( (cf_solution_location_0).c_str() );
       reader_is_open = discrete_function_reader_0.open();

       if (reader_is_open)
        { discrete_function_reader_0.read( 0, conservative_flux_coarse_ent_e0 ); }

       // flux for e_1 ...

       char location_lps_1[50];
       sprintf( location_lps_1, "_conservativeFlux_e_%d_sg_%d", 1, index_coarse_entity );
       std::string location_lps_s_1( location_lps_1 );

       std :: string cf_solution_location_1;

       // the file/place, where we saved the solutions conservative flux problems problems
       cf_solution_location_1 = path_ + "/cf_problems/" + location_lps_s_1;

       reader_is_open = false;
       // reader for data file:
       DiscreteFunctionReader discrete_function_reader_1( (cf_solution_location_1).c_str() );
       reader_is_open = discrete_function_reader_1.open();

       if (reader_is_open)
        { discrete_function_reader_1.read( 0, conservative_flux_coarse_ent_e1 ); }

       DiscreteFunctionType cflux_coarse_ent_e0_host( "Conservative Flux on coarse entity for e_0", fineDiscreteFunctionSpace_ );
       DiscreteFunctionType cflux_coarse_ent_e1_host( "Conservative Flux on coarse entity for e_1", fineDiscreteFunctionSpace_ );

       subgrid_to_hostrid_function( conservative_flux_coarse_ent_e0,
                                    conservative_flux_coarse_ent_e1,
                                    cflux_coarse_ent_e0_host,
                                    cflux_coarse_ent_e1_host );

       // flux for each neighbor entity
       DiscreteFunctionType* cflux_neighbor_ent_e0_host[ 3 ];
       DiscreteFunctionType* cflux_neighbor_ent_e1_host[ 3 ];

       std :: vector < IntersectionIteratorType > coarse_face;
       RangeType coarse_face_volume[ 3 ];

       const GridPartType &coarseGridPart = specifier_.coarseSpace().gridPart();

       int local_face_index = 0;

       IntersectionIteratorType endnit = coarseGridPart.iend( coarse_entity );
       for( IntersectionIteratorType face_it = coarseGridPart.ibegin( coarse_entity ); face_it != endnit ; ++face_it)
        {

          coarse_face.push_back( face_it );

          coarse_face_volume[ local_face_index ] = face_it->geometry().volume();

          if ( face_it->neighbor() )
            {

              EntityPointerType outside_it = face_it->outside();

              int index_coarse_neighbor_entity = coarseGridLeafIndexSet.index( *outside_it );

              // --- get subgrids and load fluxes ---

              SubGridType& sub_grid_neighbor_U_T = subgrid_list_.getSubGrid( index_coarse_neighbor_entity );

              SubGridPart subGridPart_neighbor( sub_grid_neighbor_U_T );
              SubGridDiscreteFunctionSpaceType localDiscreteFunctionSpace_neighbor( subGridPart_neighbor );

              SubGridDiscreteFunctionType conservative_flux_coarse_ent_e0_neighbor( "Conservative Flux on neighbor coarse entity for e_0", localDiscreteFunctionSpace_neighbor );
              conservative_flux_coarse_ent_e0_neighbor.clear();

              SubGridDiscreteFunctionType conservative_flux_coarse_ent_e1_neighbor( "Conservative Flux on neighbor coarse entity for e_1", localDiscreteFunctionSpace_neighbor );
              conservative_flux_coarse_ent_e1_neighbor.clear();

              // --------- load local solutions -------

              char location_lps_0_neighbor[50];
              sprintf( location_lps_0_neighbor, "_conservativeFlux_e_%d_sg_%d", 0, index_coarse_neighbor_entity );
              std::string location_lps_s_0_neighbor( location_lps_0_neighbor );

              std :: string cf_solution_location_0_neighbor;

              // the file/place, where we saved the solutions conservative flux problems problems
              cf_solution_location_0_neighbor = path_ + "/cf_problems/" + location_lps_s_0_neighbor;

              reader_is_open = false;
              // reader for data file:
              DiscreteFunctionReader discrete_function_reader_0_neighbor( (cf_solution_location_0_neighbor).c_str() );
              reader_is_open = discrete_function_reader_0_neighbor.open();

              if (reader_is_open)
               { discrete_function_reader_0_neighbor.read( 0, conservative_flux_coarse_ent_e0_neighbor ); }

              // flux for e_1 ...

              char location_lps_1_neighbor[50];
              sprintf( location_lps_1_neighbor, "_conservativeFlux_e_%d_sg_%d", 1, index_coarse_neighbor_entity );
              std::string location_lps_s_1_neighbor( location_lps_1_neighbor );

              std :: string cf_solution_location_1_neighbor;

              // the file/place, where we saved the solutions conservative flux problems problems
              cf_solution_location_1_neighbor = path_ + "/cf_problems/" + location_lps_s_1_neighbor;

              reader_is_open = false;
              // reader for data file:
              DiscreteFunctionReader discrete_function_reader_1_neighbor( (cf_solution_location_1_neighbor).c_str() );
              reader_is_open = discrete_function_reader_1_neighbor.open();

              if (reader_is_open)
               { discrete_function_reader_1_neighbor.read( 0, conservative_flux_coarse_ent_e1_neighbor ); }

              cflux_neighbor_ent_e0_host[ local_face_index ] = new DiscreteFunctionType ( "Conservative Flux on neighbor coarse entity for e_0", fineDiscreteFunctionSpace_ );

              cflux_neighbor_ent_e1_host[ local_face_index ] = new DiscreteFunctionType ( "Conservative Flux on neighbor coarse entity for e_1", fineDiscreteFunctionSpace_ );

              subgrid_to_hostrid_function( conservative_flux_coarse_ent_e0_neighbor,
                                           conservative_flux_coarse_ent_e1_neighbor,
                                           *cflux_neighbor_ent_e0_host[ local_face_index ],
                                           *cflux_neighbor_ent_e1_host[ local_face_index ] );

            }

          local_face_index += 1;
        }

       if ( local_face_index != 3 )
        { std :: cout << "Error! Implementation only for triangular mesh in 2d!" << std :: endl; abort(); }


       SubGridIteratorType sub_endit = localDiscreteFunctionSpace.end();
       for( SubGridIteratorType sub_it = localDiscreteFunctionSpace.begin(); sub_it != sub_endit; ++sub_it )
          {

             const SubGridEntityType &sub_entity = *sub_it;

             EntityPointerType host_entity_pointer = sub_grid_U_T.template getHostEntity<0>( *sub_it );
             const EntityType& host_entity = *host_entity_pointer;
	     
             EntityPointerType father_of_sub_grid_entity = host_entity_pointer;
             for (int lev = 0; lev < specifier_.getLevelDifference() ; ++lev)
               father_of_sub_grid_entity = father_of_sub_grid_entity->father();

             bool father_found = coarseGridLeafIndexSet.contains( *father_of_sub_grid_entity );
             while ( father_found == false )
               {
                 father_of_sub_grid_entity = father_of_sub_grid_entity->father();
                 father_found = coarseGridLeafIndexSet.contains( *father_of_sub_grid_entity );
               }
          
             int coarse_sub_father_index = coarseGridLeafIndexSet.index( *father_of_sub_grid_entity );
	     if ( coarse_sub_father_index != index_coarse_entity )
	      { continue; }
	      
             LocalFunctionType loc_cf_coarse_ent_e0 = cflux_coarse_ent_e0_host.localFunction( host_entity );
             LocalFunctionType loc_cf_coarse_ent_e1 = cflux_coarse_ent_e1_host.localFunction( host_entity );

             LocalFunctionType loc_msfem_coarse_part = msfem_coarse_part.localFunction( host_entity );


             IntersectionIteratorType end_it_U_T = fineDiscreteFunctionSpace_.gridPart().iend( host_entity );
             for( IntersectionIteratorType face_it_U_T = fineDiscreteFunctionSpace_.gridPart().ibegin( host_entity );
                  face_it_U_T != end_it_U_T ; ++face_it_U_T)
               {

                int relevant_face_index = -1;

                if ( is_subface( *face_it_U_T, *coarse_face[0] ) )
                  { relevant_face_index = 0; }

                if ( is_subface( *face_it_U_T, *coarse_face[1] ) )
                  { relevant_face_index = 1; }

                if ( is_subface( *face_it_U_T, *coarse_face[2] ) )
                  { relevant_face_index = 2; }

                if ( ( relevant_face_index == -1 ) || ( face_it_U_T->neighbor() == false ) )
                 { continue; }

                EntityPointerType outside_sub_it = face_it_U_T->outside();

                LocalFunctionType loc_cf_coarse_neighbor_ent_e0 = (*cflux_neighbor_ent_e0_host[ relevant_face_index ]).localFunction( host_entity );
                LocalFunctionType loc_cf_coarse_neighbor_ent_e1 = (*cflux_neighbor_ent_e1_host[ relevant_face_index ]).localFunction( host_entity );


                LocalFunctionType loc_msfem_coarse_part_neighbor = msfem_coarse_part.localFunction( *outside_sub_it );

                // evaluate the gradient of the MsfEM coarse part in the center of the coarse entity
                EntityQuadratureType coarseEntQuadrature( host_entity , 0 );
                JacobianRangeType gradient_msfem_coarse_ent(0.);
                loc_msfem_coarse_part.jacobian(coarseEntQuadrature[ 0 ], gradient_msfem_coarse_ent );

                // evaluate the gradient of the MsfEM coarse part in the center of the current neighbor of the coarse entity
                EntityQuadratureType coarseEntQuadratureNeighbor( *outside_sub_it , 0 );
                JacobianRangeType gradient_msfem_coarse_neighbor_ent(0.);
                loc_msfem_coarse_part_neighbor.jacobian(coarseEntQuadratureNeighbor[ 0 ], gradient_msfem_coarse_neighbor_ent );


                FaceQuadratureType faceQuadrature( fineDiscreteFunctionSpace_.gridPart(), *face_it_U_T, 2*fineDiscreteFunctionSpace_.order()+2 , FaceQuadratureType::INSIDE);
                 // inside macht hier keinen Unterschied, da wir formal stetige Funktionen haben und nicht die Gradienten auswerten

                 const FaceGeometryType& faceGeometry = face_it_U_T->geometry();

                 const size_t numQuadraturePoints = faceQuadrature.nop();

                 RangeType jump_integral( 0.0 );

                 RangeType check_sum( 0.0 );
                 for( size_t quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint )
                   {

                      typedef typename Intersection::LocalCoordinate LocalCoordinate;
                      const LocalCoordinate local_point = faceGeometry.local( faceQuadrature.point( quadraturePoint ) );

                      // integration factors
                      const double integrationFactor = faceGeometry.integrationElement( local_point );

                      // weight
                      const double quadratureWeight = faceQuadrature.weight(  quadraturePoint );

                      check_sum += integrationFactor * quadratureWeight;


                      RangeType value_ent_e0;
                      loc_cf_coarse_ent_e0.evaluate( faceQuadrature[ quadraturePoint ], value_ent_e0 );

                      RangeType value_neighbor_ent_e0;
                      loc_cf_coarse_neighbor_ent_e0.evaluate( faceQuadrature[ quadraturePoint ], value_neighbor_ent_e0 );

                      RangeType value_ent_e1;
                      loc_cf_coarse_ent_e1.evaluate( faceQuadrature[ quadraturePoint ], value_ent_e1 );

                      RangeType value_neighbor_ent_e1;
                      loc_cf_coarse_neighbor_ent_e1.evaluate( faceQuadrature[ quadraturePoint ], value_neighbor_ent_e1 );

                      RangeType H_E = coarse_face_volume[ relevant_face_index ];

                      // wenn man tunen moechtet... fabs() einbinden und das Vorzeichen davor aendern.

                      RangeType jump_contribution = gradient_msfem_coarse_ent[0][0] * ( /*fabs*/(value_ent_e0) + /*-fabs*/(value_neighbor_ent_e0));
                      jump_contribution +=  gradient_msfem_coarse_ent[0][1] * (/*fabs*/(value_ent_e1) + /*-fabs*/(value_neighbor_ent_e1));

                      jump[ relevant_face_index ] += H_E * integrationFactor * quadratureWeight * pow( jump_contribution, 2.0 );


                      jump_contribution = gradient_msfem_coarse_neighbor_ent[0][0] * (/*fabs*/(value_ent_e0) + /*-fabs*/(value_neighbor_ent_e0));
                      jump_contribution += gradient_msfem_coarse_neighbor_ent[0][1] * (/*fabs*/(value_ent_e1) + /*-fabs*/(value_neighbor_ent_e1));

                      jump[ relevant_face_index ] += H_E * integrationFactor * quadratureWeight * pow( jump_contribution, 2.0 );


                      JacobianRangeType coarse_jump_contribution;
                      coarse_jump_contribution[0] = gradient_msfem_coarse_ent[0] - gradient_msfem_coarse_neighbor_ent[0];

     
                      RangeType cjump(0.0);
                      cjump += fabs( value_ent_e0 * coarse_jump_contribution[0][0] + value_ent_e1 * coarse_jump_contribution[0][1] );
                      cjump += fabs( value_neighbor_ent_e0 * coarse_jump_contribution[0][0] + value_neighbor_ent_e1 * coarse_jump_contribution[0][1] );
		      
                      coarse_jump[ relevant_face_index ] += H_E * integrationFactor * quadratureWeight * pow( cjump , 2.0 );
   
#if 0
std :: cout << "index_coarse_entity = " << index_coarse_entity << std ::endl;
std :: cout << "(local_point,quadraturePoint) = (" << local_point << "," << quadraturePoint << ")" << std :: endl;

std :: cout << "value_ent_e0 = " << value_ent_e0 << std :: endl;
std :: cout << "value_neighbor_ent_e0 = " << value_neighbor_ent_e0 << std :: endl;
std :: cout << "value_ent_e1 = " << value_ent_e1 << std :: endl;
std :: cout << "value_neighbor_ent_e1 = " << value_neighbor_ent_e1 << std :: endl << std :: endl;
#endif

                   } // done loop over all quadrature points
//std :: cout << "--------------------------------------------" << std :: endl << std :: endl;

                 if ( check_sum != faceGeometry.volume())
                  { std :: cout << "Error in Face Quadrature." << std :: endl; abort(); }


               }

          }



      // std :: cout << "jump[0] = " << jump[0] << std :: endl;
      // std :: cout << "jump[1] = " << jump[1] << std :: endl;
      // std :: cout << "jump[2] = " << jump[2] << std :: endl;
      // std :: cout << std :: endl;

// laufe jetzt ueber das Subgrid der Coarse entity
// wandele die sub entities in host entities um und nutze den 'subgrid intersection iterator'
// checke ob die Intersection Teil einer coarse intersection ist
// wenn ja addiere die jeweiligen conservative fluxes ausgewertet im Quadraturpunkt auf auf dem subgrid face.


      jump_conservative_flux = ( sqrt(jump[0]) + sqrt(jump[1]) + sqrt(jump[2]) );
      jump_coarse_flux = sqrt(coarse_jump[0]) + sqrt(coarse_jump[1]) + sqrt(coarse_jump[2]);

     }
#endif

#if 0

    // \eta_T^{app}
    RangeType indicator_app_1( const EntityType &entity,
                                     DiscreteFunctionType &u_H,
                               const PeriodicDiscreteFunctionType &corrector_u_H_on_entity )
     {

       RangeType local_indicator(0.0);

       Problem::ModelProblemData model_info;
       const double delta = model_info.getDelta();
       const double epsilon_estimated = model_info.getEpsilonEstimated();

       EntityQuadratureType entityQuadrature( entity , 0 ); // 0 = polynomial order
       // the global quadrature (quadrature on the macro element T)

       const EntityGeometryType& globalEntityGeometry = entity.geometry();

       const DomainType &x_T = globalEntityGeometry.global(entityQuadrature.point(0));

       const RangeType entityVolume = entityQuadrature.weight(0) * 
               globalEntityGeometry.integrationElement(entityQuadrature.point(0));

       // int_Y A_h^{\epsilon} - A^{\epsilob}
       RangeType difference[dimension];
       for( int k = 0; k < dimension; ++k )
         difference[k] = 0.0;

       // \nabla u_H(x_T)
       LocalFunctionType u_H_local = u_H.localFunction(entity);
       JacobianRangeType gradient_u_H(0.);
       u_H_local.jacobian(entityQuadrature[ 0 ], gradient_u_H);

       // iterator over the elements of the periodic micro grid:
       PeriodicIteratorType p_endit = periodicDiscreteFunctionSpace_.end();
       for( PeriodicIteratorType p_it = periodicDiscreteFunctionSpace_.begin(); p_it != p_endit; ++p_it )
        {

          const PeriodicEntityType &micro_entity = *p_it;

          int quadOrder = 2 * PeriodicDiscreteFunctionSpaceType :: polynomialOrder + 2;

          // two quadrature formulas ( A_h^{\eps}(y)=A^{\eps}(y_s) vs. A^{\eps}(y) )
          PeriodicEntityQuadratureType one_point_quadrature( micro_entity , 0 );
          PeriodicEntityQuadratureType high_order_quadrature( micro_entity , quadOrder );

          // Q_h(u_H)(x_T,y) on the micro entity:
          PeriodicLocalFunctionType loc_Q_u_H_x_T = corrector_u_H_on_entity.localFunction(micro_entity);
          JacobianRangeType gradient_Q_u_H_x_T(0.);
          loc_Q_u_H_x_T.jacobian(one_point_quadrature[ 0 ], gradient_Q_u_H_x_T);

          // S denotes the micro grid element (i.e. 'micro_entity')
          const PeriodicEntityGeometryType& geometry_S = micro_entity.geometry();
          // y_S denotes the barycenter of the micro grid element S:
          const DomainType &y_S = geometry_S.global(one_point_quadrature.point(0));

          // to evaluate A^{\epsilon}_h:
          DomainType globalPoint_center;
          JacobianRangeType direction_of_diffusion;
          for( int k = 0; k < dimension; ++k )
            {
             globalPoint_center[ k ] = x_T[ k ] + (delta * y_S[ k ]);
             direction_of_diffusion[ 0 ][ k ] = gradient_u_H[ 0 ][ k ] + gradient_Q_u_H_x_T[ 0 ][ k ];
            }

          JacobianRangeType diffusive_flux_y_S;
          diffusion_.diffusiveFlux( globalPoint_center,
                                    direction_of_diffusion,
                                    diffusive_flux_y_S );

          double cutting_function = 1.0;
          for( int k = 0; k < dimension; ++k )
            {
              // is the current quadrature point in the relevant cell?
              // Y = [ -0.5 , 0.5 ]^dimension
              if ( fabs(y_S[ k ]) > (0.5*(epsilon_estimated/delta)) )
                { cutting_function *= 0.0; }
            }

          const size_t numQuadraturePoints = high_order_quadrature.nop();
          for( size_t microQuadraturePoint = 0; microQuadraturePoint < numQuadraturePoints; ++microQuadraturePoint )
            {

              // local (barycentric) coordinates (with respect to entity)
              const typename PeriodicEntityQuadratureType::CoordinateType &y_barycentric = high_order_quadrature.point( microQuadraturePoint );

              // current quadrature point y in global coordinates
              DomainType y = geometry_S.global( y_barycentric );

              const double weight_micro_quadrature = high_order_quadrature.weight(  microQuadraturePoint ) * geometry_S.integrationElement( y_barycentric );

              // to evaluate A^{\epsilon}:
              DomainType globalPoint;
              for( int k = 0; k < dimension; ++k )
                globalPoint[ k ] = (delta * y[ k ]) + x_T[ k ];

              // diffusive flux in y (in direction \nabla u_H(x_T) + \nabla_y Q_h^T(u_H)(y_S) )
              JacobianRangeType diffusive_flux_y;
              diffusion_.diffusiveFlux( globalPoint,
                                        direction_of_diffusion,
                                        diffusive_flux_y );

              for( int k = 0; k < dimension; ++k )
                difference[ k ] += cutting_function * weight_micro_quadrature * ( diffusive_flux_y[ 0 ][ k ] - diffusive_flux_y_S[ 0 ][ k ] );

            }

        }

       for( int k = 0; k < dimension; ++k )
          local_indicator += pow( difference[ k ], 2.0 ) * entityVolume;

       return local_indicator;

     } // end of method


    // \bar{\eta}_T^{app}
    RangeType indicator_app_2( const EntityType &entity,
                                     DiscreteFunctionType &u_H,
                               const PeriodicDiscreteFunctionType &corrector_u_H_on_entity )
     {

       RangeType local_indicator(0.0);

       Problem::ModelProblemData model_info;
       const double delta = model_info.getDelta();
       const double epsilon_estimated = model_info.getEpsilonEstimated();

       EntityQuadratureType entityQuadrature( entity , 0 ); // 0 = polynomial order
       // the global quadrature (quadrature on the macro element T)

       const EntityGeometryType& globalEntityGeometry = entity.geometry();

       const DomainType &x_T = globalEntityGeometry.global(entityQuadrature.point(0));

       const RangeType entityVolume = entityQuadrature.weight(0) * 
               globalEntityGeometry.integrationElement(entityQuadrature.point(0));

       // int_Y A_h^{\epsilon} - A^{\epsilob}
       RangeType difference[dimension];
       for( int k = 0; k < dimension; ++k )
         difference[k] = 0.0;

       // \nabla u_H(x_T)
       LocalFunctionType u_H_local = u_H.localFunction(entity);
       JacobianRangeType gradient_u_H(0.);
       u_H_local.jacobian(entityQuadrature[ 0 ], gradient_u_H);

       // iterator over the elements of the periodic micro grid:
       PeriodicIteratorType p_endit = periodicDiscreteFunctionSpace_.end();
       for( PeriodicIteratorType p_it = periodicDiscreteFunctionSpace_.begin(); p_it != p_endit; ++p_it )
        {

          const PeriodicEntityType &micro_entity = *p_it;

          int quadOrder = 2 * PeriodicDiscreteFunctionSpaceType :: polynomialOrder + 2;

          // two quadrature formulas ( A_h^{\eps}(y)=A^{\eps}(y_s) vs. A^{\eps}(y) )
          PeriodicEntityQuadratureType one_point_quadrature( micro_entity , 0 );
          PeriodicEntityQuadratureType high_order_quadrature( micro_entity , quadOrder );

          // Q_h(u_H)(x_T,y) on the micro entity:
          PeriodicLocalFunctionType loc_Q_u_H_x_T = corrector_u_H_on_entity.localFunction(micro_entity);
          JacobianRangeType gradient_Q_u_H_x_T(0.);
          loc_Q_u_H_x_T.jacobian(one_point_quadrature[ 0 ], gradient_Q_u_H_x_T);

          // S denotes the micro grid element (i.e. 'micro_entity')
          const PeriodicEntityGeometryType& geometry_S = micro_entity.geometry();
          // y_S denotes the barycenter of the micro grid element S:
          const DomainType &y_S = geometry_S.global(one_point_quadrature.point(0));

          // to evaluate A^{\epsilon}_h:
          DomainType globalPoint_center;
          JacobianRangeType direction_of_diffusion;
          for( int k = 0; k < dimension; ++k )
            {
             globalPoint_center[ k ] = x_T[ k ] + (delta * y_S[ k ]);
             direction_of_diffusion[ 0 ][ k ] = gradient_u_H[ 0 ][ k ] + gradient_Q_u_H_x_T[ 0 ][ k ];
            }

          JacobianRangeType diffusive_flux_y_S;
          diffusion_.diffusiveFlux( globalPoint_center,
                                    direction_of_diffusion,
                                    diffusive_flux_y_S );

          const size_t numQuadraturePoints = high_order_quadrature.nop();
          for( size_t microQuadraturePoint = 0; microQuadraturePoint < numQuadraturePoints; ++microQuadraturePoint )
            {

              // local (barycentric) coordinates (with respect to entity)
              const typename PeriodicEntityQuadratureType::CoordinateType &y_barycentric = high_order_quadrature.point( microQuadraturePoint );

              // current quadrature point y in global coordinates
              DomainType y = geometry_S.global( y_barycentric );

              const double weight_micro_quadrature = high_order_quadrature.weight(  microQuadraturePoint ) * geometry_S.integrationElement( y_barycentric );

              // to evaluate A^{\epsilon}:
              DomainType globalPoint;
              for( int k = 0; k < dimension; ++k )
                globalPoint[ k ] = (delta * y[ k ]) + x_T[ k ];

              // diffusive flux in y (in direction \nabla u_H(x_T) + \nabla_y Q_h^T(u_H)(y_S) )
              JacobianRangeType diffusive_flux_y;
              diffusion_.diffusiveFlux( globalPoint,
                                        direction_of_diffusion,
                                        diffusive_flux_y );

              for( int k = 0; k < dimension; ++k )
                difference[ k ] += weight_micro_quadrature * ( diffusive_flux_y[ 0 ][ k ] - diffusive_flux_y_S[ 0 ][ k ] );

            }

        }

       for( int k = 0; k < dimension; ++k )
          local_indicator += pow( difference[ k ], 2.0 ) * entityVolume;

       return local_indicator;

     } // end of method



    // (1/2)*H_E^3*\eta_E^{res}
    // (1/2) since we visit every entity twice (inner and outer)
    template< class IntersectionType >
    RangeType indicator_res_E( const IntersectionType &intersection,
                                     DiscreteFunctionType &u_H,
                               const PeriodicDiscreteFunctionType &corrector_u_H_on_inner_entity,
                               const PeriodicDiscreteFunctionType &corrector_u_H_on_outer_entity )
     {

       RangeType local_indicator(0.0);

       Problem::ModelProblemData model_info;
       const double delta = model_info.getDelta();
       const double epsilon_estimated = model_info.getEpsilonEstimated();


       EntityPointerType inner_it =  intersection.inside();
       EntityPointerType outer_it =  intersection.outside();

       const EntityType& inner_entity = *inner_it;
       const EntityType& outer_entity = *outer_it;

       FaceQuadratureType faceQuadrature( discreteFunctionSpace_.gridPart(), intersection, 0 , FaceQuadratureType::INSIDE);

       DomainType unitOuterNormal =
                intersection.unitOuterNormal(faceQuadrature.localPoint(0));

       const FaceGeometryType& faceGeometry = intersection.geometry();
       // H_E (= |E|):
       const RangeType edge_length = faceGeometry.volume();

       // jump = innerValue - outerValue
       RangeType innerValue = 0.0;
       RangeType outerValue = 0.0;

       EntityQuadratureType innerEntityQuadrature( inner_entity , 0 );
       EntityQuadratureType outerEntityQuadrature( outer_entity , 0 ); // 0 = polynomial order
       // the global quadrature (quadrature on the macro element T)

       const EntityGeometryType& globalInnerEntityGeometry = inner_entity.geometry();
       const EntityGeometryType& globalOuterEntityGeometry = outer_entity.geometry();

       const DomainType &x_T_inner = globalInnerEntityGeometry.global(innerEntityQuadrature.point(0));
       const DomainType &x_T_outer = globalOuterEntityGeometry.global(outerEntityQuadrature.point(0));

       const RangeType innerEntityVolume = innerEntityQuadrature.weight(0) * 
               globalInnerEntityGeometry.integrationElement(innerEntityQuadrature.point(0));
       const RangeType outerEntityVolume = outerEntityQuadrature.weight(0) * 
               globalOuterEntityGeometry.integrationElement(outerEntityQuadrature.point(0));

       // \nabla u_H(x_T) (on the inner element T)
       LocalFunctionType inner_u_H_local = u_H.localFunction(inner_entity);
       JacobianRangeType gradient_inner_u_H(0.);
       inner_u_H_local.jacobian(innerEntityQuadrature[ 0 ], gradient_inner_u_H);

       // \nabla u_H(x_T) (on the outer element \bar{T})
       LocalFunctionType outer_u_H_local = u_H.localFunction(outer_entity);
       JacobianRangeType gradient_outer_u_H(0.);
       outer_u_H_local.jacobian(outerEntityQuadrature[ 0 ], gradient_outer_u_H);

       // iterator over the elements of the periodic micro grid:
       PeriodicIteratorType p_endit = periodicDiscreteFunctionSpace_.end();
       for( PeriodicIteratorType p_it = periodicDiscreteFunctionSpace_.begin(); p_it != p_endit; ++p_it )
        {

          const PeriodicEntityType &micro_entity = *p_it;


          //one point quadrature formula (since A_h^{\eps}(y)=A^{\eps}(y_s) )
          PeriodicEntityQuadratureType one_point_quadrature( micro_entity , 0 );

          // Q_h(u_H)(x_T,y) on the micro entity:
          PeriodicLocalFunctionType loc_Q_u_H_x_T_inner = corrector_u_H_on_inner_entity.localFunction(micro_entity);
          JacobianRangeType gradient_Q_u_H_x_T_inner(0.);
          loc_Q_u_H_x_T_inner.jacobian(one_point_quadrature[ 0 ], gradient_Q_u_H_x_T_inner);

          // Q_h(u_H)(x_{bar{T}},y) on the micro entity:
          PeriodicLocalFunctionType loc_Q_u_H_x_T_outer = corrector_u_H_on_outer_entity.localFunction(micro_entity);
          JacobianRangeType gradient_Q_u_H_x_T_outer(0.);
          loc_Q_u_H_x_T_outer.jacobian(one_point_quadrature[ 0 ], gradient_Q_u_H_x_T_outer);

          // S denotes the micro grid element (i.e. 'micro_entity')
          const PeriodicEntityGeometryType& geometry_S = micro_entity.geometry();

          // local (barycentric) coordinates (with respect to entity)
          const typename PeriodicEntityQuadratureType::CoordinateType &y_S_barycentric = one_point_quadrature.point(0);

          // y_S denotes the barycenter of the micro grid element S:
          const DomainType &y_S = geometry_S.global(y_S_barycentric);

          // to evaluate A^{\epsilon}_h:
          DomainType globalPoint_inner;
          JacobianRangeType direction_of_diffusion_inner;
          for( int k = 0; k < dimension; ++k )
           {
             globalPoint_inner[ k ] = x_T_inner[ k ] + (delta * y_S[ k ]);
             direction_of_diffusion_inner[ 0 ][ k ] = gradient_inner_u_H[ 0 ][ k ] + gradient_Q_u_H_x_T_inner[ 0 ][ k ];
           }

          // to evaluate A^{\epsilon}_h:
          DomainType globalPoint_outer;
          JacobianRangeType direction_of_diffusion_outer;
          for( int k = 0; k < dimension; ++k )
            {
             globalPoint_outer[ k ] = x_T_outer[ k ] + (delta * y_S[ k ]);
             direction_of_diffusion_outer[ 0 ][ k ] = gradient_outer_u_H[ 0 ][ k ] + gradient_Q_u_H_x_T_outer[ 0 ][ k ];
            }


          JacobianRangeType diffusive_flux_y_S_inner;
          diffusion_.diffusiveFlux( globalPoint_inner,
                                    direction_of_diffusion_inner,
                                    diffusive_flux_y_S_inner );

          JacobianRangeType diffusive_flux_y_S_outer;
          diffusion_.diffusiveFlux( globalPoint_outer,
                                    direction_of_diffusion_outer,
                                    diffusive_flux_y_S_outer );

          double cutting_function = 1.0;
          for( int k = 0; k < dimension; ++k )
            {
              // is the current quadrature point in the relevant cell?
              // Y = [ -0.5 , 0.5 ]^dimension
              if ( fabs(y_S[ k ]) > (0.5*(epsilon_estimated/delta)) )
                { cutting_function *= 0.0; }
            }

          const double weight_micro_quadrature = one_point_quadrature.weight(0) * geometry_S.integrationElement( y_S_barycentric );

          for( int k = 0; k < dimension; ++k )
            {
             innerValue += cutting_function * weight_micro_quadrature * diffusive_flux_y_S_inner[ 0 ][ k ] * unitOuterNormal[ k ];
             outerValue += cutting_function * weight_micro_quadrature * diffusive_flux_y_S_outer[ 0 ][ k ] * unitOuterNormal[ k ];
            }

        }

       local_indicator += pow( innerValue - outerValue, 2.0 );

       // in 2D (and this is what we assume) it holds |E|*H_E^3 = H_E^4 = |E|^4
       return ( pow(edge_length,4.0 ) * (1.0/2.0) * local_indicator );

     } // end of method


    // NOTE: This method ONLY works for uniformly refined micro-grid and GRIDDIM=2!
    // (even though it might be extended easily to other cases, here we use the cheapest strategy, which just works for the case described above!)
    // here, uniform also means that we have the same number of elements in every direction!)
    // \bar{\eta}_T^{res}
    RangeType indicator_res_T( const EntityType &entity,
                                     DiscreteFunctionType &u_H,
                               const PeriodicDiscreteFunctionType &corrector_u_H_on_entity )
     {


       RangeType local_indicator(0.0);

       Problem::ModelProblemData model_info;
       const double delta = model_info.getDelta();
       const double epsilon_estimated = model_info.getEpsilonEstimated();

       EntityQuadratureType entityQuadrature( entity , 0 ); // 0 = polynomial order
       // the global quadrature (quadrature on the macro element T)

       const EntityGeometryType& globalEntityGeometry = entity.geometry();

       const DomainType &x_T = globalEntityGeometry.global(entityQuadrature.point(0));

       const RangeType entityVolume = entityQuadrature.weight(0) * 
               globalEntityGeometry.integrationElement(entityQuadrature.point(0));

       // \nabla u_H(x_T)
       LocalFunctionType u_H_local = u_H.localFunction(entity);
       JacobianRangeType gradient_u_H(0.);
       u_H_local.jacobian(entityQuadrature[ 0 ], gradient_u_H);


       if ( dimension != 2 )
        {
          std :: cout << "The error indicator 'indicator_res_T' is not implemented for dimension!=2 and only works for uniformly refined micro-grids!" << std :: endl;
        }

       // edge length of a boundary face
       RangeType ref_edge_length = 1.0;

       GridPartType auxGridPart = auxiliaryDiscreteFunctionSpace_.gridPart();

       IteratorType micro_it = auxiliaryDiscreteFunctionSpace_.begin();

       // we just need one element to determine all the properties (due to uniform refinement)!
       IntersectionIteratorType endnit = auxGridPart.iend(*micro_it);
       for(IntersectionIteratorType nit = auxGridPart.ibegin(*micro_it); nit != endnit ; ++nit)
         {
           FaceQuadratureType faceQuadrature( auxGridPart, *nit, 1 , FaceQuadratureType::INSIDE);
           const FaceGeometryType& faceGeometry = nit->geometry();

           if ( ref_edge_length > faceGeometry.volume() )
            { ref_edge_length = faceGeometry.volume(); }

         }

       // number of boundary faces per cube-edge:
       int num_boundary_faces_per_direction = int( (1 / ref_edge_length) + 0.2 );
       // (+0.2 to avoid rounding errors)

       // generalized jump up/down
       RangeType jump_up_down[ num_boundary_faces_per_direction ];

       // generalized jump left/right
       RangeType jump_left_right[ num_boundary_faces_per_direction ];


       for( int id = 0; id < num_boundary_faces_per_direction; ++id )
         {
          jump_up_down[ id ] = 0.0;
          jump_left_right[ id ] = 0.0;
         }


       // procedure for computing the generalized jumps:
       // 1. compute the center of the current face E: (x_E,y_E)
       // 2. one and only one of these values fulfiles abs()=1/2
       //    (i.e. we either have 'abs(x_E)=1/2' or 'abs(y_E)=1/2')
       //    without loss of generality we assume 'abs(y_E)=1/2', then
       //    we are in the setting of 'jump_up_down'
       // 3. Any 'jump_up_down' recieves a unique ID by
       //     jump up down over (x,-1/2) and (x,1/2) = jump_up_down[ (num_boundary_faces_per_direction/2) + (( x - (edge_length/2) ) / edge_length) ]
       //     (num_boundary_faces_per_direction/2) + (( x - (edge_length/2) ) / edge_length) is the corresponding ID
       // 4. the situation is identical for 'abs(x_E)=1/2'.


       // iterator over the elements of the periodic micro grid:
       IteratorType p_endit = auxiliaryDiscreteFunctionSpace_.end();
       for( IteratorType p_it = auxiliaryDiscreteFunctionSpace_.begin(); p_it != p_endit; ++p_it )
        {

          // --------- the 'inner entity' (micro grid) ------------

          const EntityType &micro_entity = *p_it;

          // one point quadrature formula ( A_h^{\eps}(y)=A^{\eps}(y_s) )
          EntityQuadratureType one_point_quadrature( micro_entity , 0 );

          // Q_h(u_H)(x_T,y) on the micro entity:
          PeriodicLocalFunctionType loc_Q_u_H_x_T = corrector_u_H_on_entity.localFunction(micro_entity);
          JacobianRangeType gradient_Q_u_H_x_T(0.);
          loc_Q_u_H_x_T.jacobian(one_point_quadrature[ 0 ], gradient_Q_u_H_x_T);

          // S denotes the micro grid element (i.e. 'micro_entity')
          const EntityGeometryType& geometry_S = micro_entity.geometry();
          // y_S denotes the barycenter of the micro grid element S:
          const DomainType &y_S = geometry_S.global(one_point_quadrature.point(0));

          // to evaluate A^{\epsilon}_h (in center of current inner entity):
          DomainType globalPoint;
          JacobianRangeType direction_of_diffusion;
          for( int k = 0; k < dimension; ++k )
           {
            direction_of_diffusion[ 0 ][ k ] = gradient_u_H[ 0 ][ k ] + gradient_Q_u_H_x_T[ 0 ][ k ];
            globalPoint[ k ] = x_T[ k ] + (delta * y_S[ k ]);
           }

          JacobianRangeType diffusive_flux_y_S;
          diffusion_.diffusiveFlux( globalPoint,
                                    direction_of_diffusion,
                                    diffusive_flux_y_S );

          // ----------------------------------------------------

          IntersectionIteratorType p_endnit = auxGridPart.iend(micro_entity);
          for(IntersectionIteratorType p_nit = auxGridPart.ibegin(micro_entity); p_nit != p_endnit ; ++p_nit)
            {

               // Note: we are on the zero-centered unit cube! (That's why everything works!)

               FaceQuadratureType faceQuadrature( auxGridPart, *p_nit, 0, FaceQuadratureType::INSIDE );
               const FaceGeometryType& faceGeometry = p_nit->geometry();

               const RangeType edge_length = faceGeometry.volume();

               DomainType unitOuterNormal =
                       p_nit->unitOuterNormal(faceQuadrature.localPoint(0));

               //if there is a neighbor entity (the normal gradient jumps)
               if ( p_nit->neighbor() )
                {

                  // --------- the 'outer entity' (micro grid) ------------
                  //                ( neighbor entity )

                  EntityPointerType outer_p_it = p_nit->outside();
                  const EntityType &outer_micro_entity = *outer_p_it;

                  // one point quadrature formula ( A_h^{\eps}(y)=A^{\eps}(y_s) )
                  EntityQuadratureType outer_one_point_quadrature( outer_micro_entity , 0 );

                  // Q_h(u_H)(x_T,y) on the neighbor entity:
                  PeriodicLocalFunctionType outer_loc_Q_u_H_x_T = corrector_u_H_on_entity.localFunction(outer_micro_entity);
                  JacobianRangeType gradient_outer_Q_u_H_x_T(0.);
                  outer_loc_Q_u_H_x_T.jacobian(outer_one_point_quadrature[ 0 ], gradient_outer_Q_u_H_x_T);

                  // S denotes the micro grid element (i.e. 'micro_entity')
                  const EntityGeometryType& outer_geometry_S = outer_micro_entity.geometry();
                  // outer_y_S denotes the barycenter of the neighbor micro grid element of S:
                  const DomainType &outer_y_S = outer_geometry_S.global(outer_one_point_quadrature.point(0));

                  // to evaluate A^{\epsilon}_h (in center of current outer entity):
                  DomainType outer_globalPoint;
                  JacobianRangeType direction_of_diffusion_outside;
                  for( int k = 0; k < dimension; ++k )
                    {
                     direction_of_diffusion_outside[ 0 ][ k ] = gradient_u_H[ 0 ][ k ] + gradient_outer_Q_u_H_x_T[ 0 ][ k ];
                     outer_globalPoint[ k ] = x_T[ k ] + (delta * outer_y_S[ k ]);
                    }

                  JacobianRangeType diffusive_flux_y_S_outside;
                  diffusion_.diffusiveFlux( outer_globalPoint,
                                            direction_of_diffusion_outside,
                                            diffusive_flux_y_S_outside );

                  // ----------------------------------------------------

                  RangeType aux_value = 0.0;

                  for( int k = 0; k < dimension; ++k )
                    aux_value += ( diffusive_flux_y_S_outside[ 0 ][ k ] - diffusive_flux_y_S[ 0 ][ k ] ) * unitOuterNormal[ k ] ;

                  local_indicator += pow( edge_length , 4.0 ) * pow( aux_value , 2.0 );


                }
               else //if there is no neighbor entity, the face is a boundary face and we use the generalized gradient jumps.
                {

                  // Remember:
                  // procedure for computing the generalized jumps:
                  // 1. compute the center of the current face E: (x_E,y_E)
                  // 2. one and only one of these values fulfiles abs()=1/2
                  //    (i.e. we either have 'abs(x_E)=1/2' or 'abs(y_E)=1/2')
                  //    without loss of generality we assume 'abs(y_E)=1/2', then
                  //    we are in the setting of 'jump_up_down'
                  // 3. Any 'jump_up_down' recieves a unique ID by
                  //     jump up down over (x,-1/2) and (x,1/2) = ... jump_up_down[ (( x - (edge_length/2) ) / edge_length) ]
                  //     ... (( x - (edge_length/2) ) / edge_length) is the corresponding ID
                  // 4. the situation is identical for 'abs(x_E)=1/2'.

                  const DomainType &edge_center = geometry_S.global(faceQuadrature.point(0));

                  if ( fabs(edge_center[0]) == 0.5 )
                   {
                    // + 0.2 to avoid rounding errors!
                    int id = int( ( num_boundary_faces_per_direction / 2 ) + (( edge_center[1] - (edge_length / 2.0) ) / edge_length) + 0.2 );

                    // unit outer normal creates the correct sign!
                    for( int k = 0; k < dimension; ++k )
                       jump_up_down[ id ] += pow( edge_length , 2.0 ) * ( diffusive_flux_y_S[ 0 ][ k ] * unitOuterNormal[ k ] );

                    //std :: cout << "edge_center = " << edge_center << std :: endl;
                    //std :: cout << "edge_length = " << edge_length << std :: endl;
                    //std :: cout << "unitOuterNormal = " << unitOuterNormal << std :: endl;
                    //std :: cout << "diffusive_flux_y_S = " << diffusive_flux_y_S[ 0 ] << std :: endl;
                    //std :: cout << "jump_up_down id = " << id << std :: endl;
                    //std :: cout << "jump_up_down[" << id << "] = " << jump_up_down[ id ] << std :: endl << std :: endl;

                   }
                  else
                   {

                    // + 0.2 to avoid rounding errors!
                    int id = int( ( num_boundary_faces_per_direction / 2 ) + (( edge_center[0] - (edge_length / 2.0) ) / edge_length) + 0.2 );

                    // unit outer normal creates the correct sign!
                    for( int k = 0; k < dimension; ++k )
                       jump_left_right[ id ] += pow( edge_length , 2.0 ) * ( diffusive_flux_y_S[ 0 ][ k ] * unitOuterNormal[ k ] );

                    //std :: cout << "edge_center = " << edge_center << std :: endl;
                    //std :: cout << "edge_length = " << edge_length << std :: endl;
                    //std :: cout << "unitOuterNormal = " << unitOuterNormal << std :: endl;
                    //std :: cout << "diffusive_flux_y_S = " << diffusive_flux_y_S[ 0 ] << std :: endl;
                    //std :: cout << "jump_left_right id = " << id << std :: endl;
                    //std :: cout << "jump_left_right[" << id << "] = " << jump_left_right[ id ] << std :: endl << std :: endl;

                   }

                }

            }


        }

       for( int id = 0; id < num_boundary_faces_per_direction; ++id )
         {
          local_indicator += ( pow( jump_up_down[ id ] , 2.0 ) + pow( jump_left_right[ id ] , 2.0 ) );
         }

       local_indicator *= entityVolume;

       return local_indicator;

     } // end of method


//only for TFR:

    // NOTE: This method ONLY works for uniformly refined micro-grid and GRIDDIM=2!
    // (even though it might be extended easily to other cases, here we use the cheapest strategy, which just works for the case described above!)
    // here, uniform also means that we have the same number of elements in every direction!)
    // \eta_T^{tfr}
    RangeType indicator_tfr_1( const EntityType &entity,
                                     DiscreteFunctionType &u_H,
                               const PeriodicDiscreteFunctionType &corrector_u_H_on_entity )
     {

       RangeType local_indicator(0.0);

       Problem::ModelProblemData model_info;
       const double delta = model_info.getDelta();
       const double epsilon_estimated = model_info.getEpsilonEstimated();

       EntityQuadratureType entityQuadrature( entity , 0 ); // 0 = polynomial order
       // the global quadrature (quadrature on the macro element T)

       const EntityGeometryType& globalEntityGeometry = entity.geometry();

       const DomainType &x_T = globalEntityGeometry.global(entityQuadrature.point(0));

       const RangeType entityVolume = entityQuadrature.weight(0) * 
               globalEntityGeometry.integrationElement(entityQuadrature.point(0));

       // \nabla u_H(x_T)
       LocalFunctionType u_H_local = u_H.localFunction(entity);
       JacobianRangeType gradient_u_H(0.);
       u_H_local.jacobian(entityQuadrature[ 0 ], gradient_u_H);


       if ( dimension != 2 )
        {
          std :: cout << "The error indicator 'indicator_tfr_1' is not implemented for dimension!=2 and only works for uniformly refined micro-grids!" << std :: endl;
        }


       // edge length of a boundary face
       RangeType ref_edge_length = 1.0;

       GridPartType auxGridPart = auxiliaryDiscreteFunctionSpace_.gridPart();

       IteratorType micro_it = auxiliaryDiscreteFunctionSpace_.begin();

       // we just need one element to determine all the properties (due to uniform refinement)!
       IntersectionIteratorType endnit = auxGridPart.iend(*micro_it);
       for(IntersectionIteratorType nit = auxGridPart.ibegin(*micro_it); nit != endnit ; ++nit)
         {
           FaceQuadratureType faceQuadrature( auxGridPart, *nit, 1 , FaceQuadratureType::INSIDE);
           const FaceGeometryType& faceGeometry = nit->geometry();

           if ( ref_edge_length > faceGeometry.volume() )
            { ref_edge_length = faceGeometry.volume(); }

         }

       // number of boundary faces per (\epsilon/\delta-scaled) cube edge:
       int num_boundary_faces_per_direction = int( ((epsilon_estimated/delta) / ref_edge_length) + 0.2 );
       // (+0.2 to avoid rounding errors)

       // generalized jump up/down
       RangeType jump_up_down[ num_boundary_faces_per_direction ];

       // generalized jump left/right
       RangeType jump_left_right[ num_boundary_faces_per_direction ];

       for( int id = 0; id < num_boundary_faces_per_direction; ++id )
         {
          jump_up_down[ id ] = 0.0;
          jump_left_right[ id ] = 0.0;
         }

       // did you find a boundary edge of the \eps-\delta-cube? (you must find it!!)
       bool eps_delta_boundary_edge_found = false;

       // procedure for computing the generalized jumps:
       // 1. compute the center of the current face E: (x_E,y_E)
       // 2. one and only one of these values fulfiles abs()=1/2
       //    (i.e. we either have 'abs(x_E)=1/2' or 'abs(y_E)=1/2')
       //    without loss of generality we assume 'abs(y_E)=1/2', then
       //    we are in the setting of 'jump_up_down'
       // 3. Any 'jump_up_down' recieves a unique ID by
       //     jump up down over (x,-1/2) and (x,1/2) = jump_up_down[ (num_boundary_faces_per_direction/2) + (( x - (edge_length/2) ) / edge_length) ]
       //     (num_boundary_faces_per_direction/2) + (( x - (edge_length/2) ) / edge_length) is the corresponding ID
       // 4. the situation is identical for 'abs(x_E)=1/2'.


       // iterator over the elements of the periodic micro grid:
       IteratorType p_endit = auxiliaryDiscreteFunctionSpace_.end();
       for( IteratorType p_it = auxiliaryDiscreteFunctionSpace_.begin(); p_it != p_endit; ++p_it )
        {

          // --------- the 'inner entity' (micro grid) ------------

          const EntityType &micro_entity = *p_it;

          // one point quadrature formula ( A_h^{\eps}(y)=A^{\eps}(y_s) )
          EntityQuadratureType one_point_quadrature( micro_entity , 0 );

          // Q_h(u_H)(x_T,y) on the micro entity:
          PeriodicLocalFunctionType loc_Q_u_H_x_T = corrector_u_H_on_entity.localFunction(micro_entity);
          JacobianRangeType gradient_Q_u_H_x_T(0.);
          loc_Q_u_H_x_T.jacobian(one_point_quadrature[ 0 ], gradient_Q_u_H_x_T);

          // S denotes the micro grid element (i.e. 'micro_entity')
          const EntityGeometryType& geometry_S = micro_entity.geometry();
          // y_S denotes the barycenter of the micro grid element S:
          const DomainType &y_S = geometry_S.global(one_point_quadrature.point(0));

          // to evaluate A^{\epsilon}_h (in center of current inner entity):
          DomainType globalPoint;
          JacobianRangeType direction_of_diffusion;
          for( int k = 0; k < dimension; ++k )
            {
             direction_of_diffusion[ 0 ][ k ] = gradient_u_H[ 0 ][ k ] + gradient_Q_u_H_x_T[ 0 ][ k ];
             globalPoint[ k ] = ( delta * y_S[ k ] ) + x_T[ k ];
            }

          JacobianRangeType diffusive_flux_y_S;
          diffusion_.diffusiveFlux( globalPoint,
                                    direction_of_diffusion,
                                    diffusive_flux_y_S );

          if ( (fabs(y_S[1]) <= ((0.5*(epsilon_estimated/delta)))) &&
               (fabs(y_S[0]) <= ((0.5*(epsilon_estimated/delta))))    )
          {
          // ----------------------------------------------------

          IntersectionIteratorType p_endnit = auxGridPart.iend(micro_entity);
          for(IntersectionIteratorType p_nit = auxGridPart.ibegin(micro_entity); p_nit != p_endnit ; ++p_nit)
            {

               // Note: we are on the zero-centered unit cube! (That's why everything works!)

               FaceQuadratureType faceQuadrature( auxGridPart, *p_nit, 0, FaceQuadratureType::INSIDE );
               const FaceGeometryType& faceGeometry = p_nit->geometry();

               const RangeType edge_length = faceGeometry.volume();

               DomainType unitOuterNormal =
                       p_nit->unitOuterNormal(faceQuadrature.localPoint(0));

               const DomainType &edge_center = geometry_S.global(faceQuadrature.point(0));

               if ( (fabs(edge_center[0]) == ((0.5*(epsilon_estimated/delta)))) || 
                    (fabs(edge_center[1]) == ((0.5*(epsilon_estimated/delta)))) )
                {

                  //we use the generalized gradient jumps.

                  if ( (fabs(edge_center[0]) == ((0.5*(epsilon_estimated/delta)))) &&
                       (fabs(edge_center[1]) < ((0.5*(epsilon_estimated/delta))))   )
                   {
                    // + 0.2 to avoid rounding errors!
                    int id = int( ( num_boundary_faces_per_direction / 2 ) + (( edge_center[1] - (edge_length / 2.0) ) / edge_length) + 0.2 );

                    // unit outer normal creates the correct sign!
                    for( int k = 0; k < dimension; ++k )
                       jump_up_down[ id ] += sqrt(edge_length) * ( diffusive_flux_y_S[ 0 ][ k ] * unitOuterNormal[ k ] );

                    //std :: cout << "edge_center = " << edge_center << std :: endl;
                    //std :: cout << "edge_length = " << edge_length << std :: endl;
                    //std :: cout << "unitOuterNormal = " << unitOuterNormal << std :: endl;
                    //std :: cout << "diffusive_flux_y_S = " << diffusive_flux_y_S[ 0 ] << std :: endl;
                    //std :: cout << "jump_up_down id = " << id << std :: endl;
                    //std :: cout << "jump_up_down[" << id << "] = " << jump_up_down[ id ] << std :: endl << std :: endl;

                    eps_delta_boundary_edge_found = true;
                    //std :: cout << "Found eps/delta boundary edge with center = " << edge_center << std :: endl;

                   }

                  if ( (fabs(edge_center[1]) == ((0.5*(epsilon_estimated/delta)))) &&
                       (fabs(edge_center[0]) < ((0.5*(epsilon_estimated/delta))))   )
                   {

                    // + 0.2 to avoid rounding errors!
                    int id = int( ( num_boundary_faces_per_direction / 2 ) + (( edge_center[0] - (edge_length / 2.0) ) / edge_length) + 0.2 );

                    // unit outer normal creates the correct sign!
                    for( int k = 0; k < dimension; ++k )
                       jump_left_right[ id ] += sqrt(edge_length) * ( diffusive_flux_y_S[ 0 ][ k ] * unitOuterNormal[ k ] );

                    //std :: cout << "edge_center = " << edge_center << std :: endl;
                    //std :: cout << "edge_length = " << edge_length << std :: endl;
                    //std :: cout << "unitOuterNormal = " << unitOuterNormal << std :: endl;
                    //std :: cout << "diffusive_flux_y_S = " << diffusive_flux_y_S[ 0 ] << std :: endl;
                    //std :: cout << "jump_left_right id = " << id << std :: endl;
                    //std :: cout << "jump_left_right[" << id << "] = " << jump_left_right[ id ] << std :: endl << std :: endl;

                    eps_delta_boundary_edge_found = true;
                    //std :: cout << "Found eps/delta boundary edge with center = " << edge_center << std :: endl;

                   }

                }

               if ( ( fabs(edge_center[0]) < (0.5*(epsilon_estimated/delta)) ) && 
                    ( fabs(edge_center[1]) < (0.5*(epsilon_estimated/delta)) )  )
                {
                   // in this situation, p_nit always has an outside neighbor entity
                   // (use the normal gradient jumps)


                   // --------- the 'outer entity' (micro grid) ------------
                   //                ( neighbor entity )

                   EntityPointerType outer_p_it = p_nit->outside();
                   const EntityType &outer_micro_entity = *outer_p_it;

                   // one point quadrature formula ( A_h^{\eps}(y)=A^{\eps}(y_s) )
                   EntityQuadratureType outer_one_point_quadrature( outer_micro_entity , 0 );

                   // Q_h(u_H)(x_T,y) on the neighbor entity:
                   PeriodicLocalFunctionType outer_loc_Q_u_H_x_T = corrector_u_H_on_entity.localFunction(outer_micro_entity);
                   JacobianRangeType gradient_outer_Q_u_H_x_T(0.);
                   outer_loc_Q_u_H_x_T.jacobian(outer_one_point_quadrature[ 0 ], gradient_outer_Q_u_H_x_T);

                   // S denotes the micro grid element (i.e. 'micro_entity')
                   const EntityGeometryType& outer_geometry_S = outer_micro_entity.geometry();
                   // outer_y_S denotes the barycenter of the neighbor micro grid element of S:
                   const DomainType &outer_y_S = outer_geometry_S.global(outer_one_point_quadrature.point(0));

                   // to evaluate A^{\epsilon}_h (in center of current outer entity):
                   DomainType outer_globalPoint;
                   JacobianRangeType direction_of_diffusion_outside;
                   for( int k = 0; k < dimension; ++k )
                     {
                      outer_globalPoint[ k ] = x_T[ k ] + (delta * outer_y_S[ k ]);
                      direction_of_diffusion_outside[ 0 ][ k ] = gradient_u_H[ 0 ][ k ] +  gradient_outer_Q_u_H_x_T[ 0 ][ k ];
                     }

                   JacobianRangeType diffusive_flux_y_S_outside;
                   diffusion_.diffusiveFlux( outer_globalPoint,
                                             direction_of_diffusion_outside,
                                             diffusive_flux_y_S_outside );

                   // ----------------------------------------------------

                   RangeType aux_value = 0.0;

                   for( int k = 0; k < dimension; ++k )
                     aux_value += ( diffusive_flux_y_S_outside[ 0 ][ k ] - diffusive_flux_y_S[ 0 ][ k ] ) * unitOuterNormal[ k ] ;

                   local_indicator += edge_length * pow( aux_value , 2.0 );
                }

            } // end Intersection iterator
        } // endif: "y_S \in \eps/\delta Y"

        }

       for( int id = 0; id < num_boundary_faces_per_direction; ++id )
         {
          local_indicator += ( pow( jump_up_down[ id ] , 2.0 ) + pow( jump_left_right[ id ] , 2.0 ) );
         }

       if ( eps_delta_boundary_edge_found == false )
        {
           std :: cout << "Error! Make sure that the restriction of the Y-triangulation on 'eps/delta Y' is a complete periodic triangulation of 'eps/delta Y' on its own. (for instance: delta = 2 epsilon should work)" << std :: endl; 
           std :: abort();
        }

       local_indicator *= entityVolume;

       return local_indicator;

     } // end of method



#if 1

    // NOTE: This method ONLY works for uniformly refined micro-grid and GRIDDIM=2!
    // (even though it might be extended easily to other cases, here we use the cheapest strategy, which just works for the case described above!)
    // here, uniform also means that we have the same number of elements in every direction!)
    // This indicator does not exist in theory. It is is additionally multiplied with h^3 to fit the other orders of convergence. In this setting it is more suitable for a comparison to capture the effect of boundary jumps in the case of a wrong boundary condition
    RangeType indicator_effective_tfr( const EntityType &entity,
                                       DiscreteFunctionType &u_H,
                                       const PeriodicDiscreteFunctionType &corrector_u_H_on_entity )
     {

       RangeType local_indicator(0.0);

       Problem::ModelProblemData model_info;
       const double delta = model_info.getDelta();
       const double epsilon_estimated = model_info.getEpsilonEstimated();

       EntityQuadratureType entityQuadrature( entity , 0 ); // 0 = polynomial order
       // the global quadrature (quadrature on the macro element T)

       const EntityGeometryType& globalEntityGeometry = entity.geometry();

       const DomainType &x_T = globalEntityGeometry.global(entityQuadrature.point(0));

       const RangeType entityVolume = entityQuadrature.weight(0) * 
               globalEntityGeometry.integrationElement(entityQuadrature.point(0));

       // \nabla u_H(x_T)
       LocalFunctionType u_H_local = u_H.localFunction(entity);
       JacobianRangeType gradient_u_H(0.);
       u_H_local.jacobian(entityQuadrature[ 0 ], gradient_u_H);


       if ( dimension != 2 )
        {
          std :: cout << "The error indicator 'indicator_tfr_1' is not implemented for dimension!=2 and only works for uniformly refined micro-grids!" << std :: endl;
        }


       // edge length of a boundary face
       RangeType ref_edge_length = 1.0;

       GridPartType auxGridPart = auxiliaryDiscreteFunctionSpace_.gridPart();

       IteratorType micro_it = auxiliaryDiscreteFunctionSpace_.begin();

       // we just need one element to determine all the properties (due to uniform refinement)!
       IntersectionIteratorType endnit = auxGridPart.iend(*micro_it);
       for(IntersectionIteratorType nit = auxGridPart.ibegin(*micro_it); nit != endnit ; ++nit)
         {
           FaceQuadratureType faceQuadrature( auxGridPart, *nit, 1 , FaceQuadratureType::INSIDE);
           const FaceGeometryType& faceGeometry = nit->geometry();

           if ( ref_edge_length > faceGeometry.volume() )
            { ref_edge_length = faceGeometry.volume(); }

         }

       // number of boundary faces per (\epsilon/\delta-scaled) cube edge:
       int num_boundary_faces_per_direction = int( ((epsilon_estimated/delta) / ref_edge_length) + 0.2 );
       // (+0.2 to avoid rounding errors)

       // generalized jump up/down
       RangeType jump_up_down[ num_boundary_faces_per_direction ];

       // generalized jump left/right
       RangeType jump_left_right[ num_boundary_faces_per_direction ];

       for( int id = 0; id < num_boundary_faces_per_direction; ++id )
         {
          jump_up_down[ id ] = 0.0;
          jump_left_right[ id ] = 0.0;
         }

       // did you find a boundary edge of the \eps-\delta-cube? (you must find it!!)
       bool eps_delta_boundary_edge_found = false;

       // procedure for computing the generalized jumps:
       // 1. compute the center of the current face E: (x_E,y_E)
       // 2. one and only one of these values fulfiles abs()=1/2
       //    (i.e. we either have 'abs(x_E)=1/2' or 'abs(y_E)=1/2')
       //    without loss of generality we assume 'abs(y_E)=1/2', then
       //    we are in the setting of 'jump_up_down'
       // 3. Any 'jump_up_down' recieves a unique ID by
       //     jump up down over (x,-1/2) and (x,1/2) = jump_up_down[ (num_boundary_faces_per_direction/2) + (( x - (edge_length/2) ) / edge_length) ]
       //     (num_boundary_faces_per_direction/2) + (( x - (edge_length/2) ) / edge_length) is the corresponding ID
       // 4. the situation is identical for 'abs(x_E)=1/2'.


       // iterator over the elements of the periodic micro grid:
       IteratorType p_endit = auxiliaryDiscreteFunctionSpace_.end();
       for( IteratorType p_it = auxiliaryDiscreteFunctionSpace_.begin(); p_it != p_endit; ++p_it )
        {

          // --------- the 'inner entity' (micro grid) ------------

          const EntityType &micro_entity = *p_it;

          // one point quadrature formula ( A_h^{\eps}(y)=A^{\eps}(y_s) )
          EntityQuadratureType one_point_quadrature( micro_entity , 0 );

          // Q_h(u_H)(x_T,y) on the micro entity:
          PeriodicLocalFunctionType loc_Q_u_H_x_T = corrector_u_H_on_entity.localFunction(micro_entity);
          JacobianRangeType gradient_Q_u_H_x_T(0.);
          loc_Q_u_H_x_T.jacobian(one_point_quadrature[ 0 ], gradient_Q_u_H_x_T);

          // S denotes the micro grid element (i.e. 'micro_entity')
          const EntityGeometryType& geometry_S = micro_entity.geometry();
          // y_S denotes the barycenter of the micro grid element S:
          const DomainType &y_S = geometry_S.global(one_point_quadrature.point(0));

          // to evaluate A^{\epsilon}_h (in center of current inner entity):
          DomainType globalPoint;
          JacobianRangeType direction_of_diffusion;
          for( int k = 0; k < dimension; ++k )
           {
            direction_of_diffusion[ 0 ][ k ] = gradient_u_H[ 0 ][ k ] + gradient_Q_u_H_x_T[ 0 ][ k ];
            globalPoint[ k ] = x_T[ k ] + (delta * y_S[ k ]);
           }


          JacobianRangeType diffusive_flux_y_S;
          diffusion_.diffusiveFlux( globalPoint,
                                    direction_of_diffusion,
                                    diffusive_flux_y_S );

          if ( (fabs(y_S[1]) <= ((0.5*(epsilon_estimated/delta)))) &&
               (fabs(y_S[0]) <= ((0.5*(epsilon_estimated/delta))))    )
          {
          // ----------------------------------------------------

          IntersectionIteratorType p_endnit = auxGridPart.iend(micro_entity);
          for(IntersectionIteratorType p_nit = auxGridPart.ibegin(micro_entity); p_nit != p_endnit ; ++p_nit)
            {

               // Note: we are on the zero-centered unit cube! (That's why everything works!)

               FaceQuadratureType faceQuadrature( auxGridPart, *p_nit, 0, FaceQuadratureType::INSIDE );
               const FaceGeometryType& faceGeometry = p_nit->geometry();

               const RangeType edge_length = faceGeometry.volume();

               DomainType unitOuterNormal =
                       p_nit->unitOuterNormal(faceQuadrature.localPoint(0));

               const DomainType &edge_center = geometry_S.global(faceQuadrature.point(0));

               if ( (fabs(edge_center[0]) == ((0.5*(epsilon_estimated/delta)))) || 
                    (fabs(edge_center[1]) == ((0.5*(epsilon_estimated/delta)))) )
                {

                  //we use the generalized gradient jumps.

                  if ( (fabs(edge_center[0]) == ((0.5*(epsilon_estimated/delta)))) &&
                       (fabs(edge_center[1]) < ((0.5*(epsilon_estimated/delta))))   )
                   {
                    // + 0.2 to avoid rounding errors!
                    int id = int( ( num_boundary_faces_per_direction / 2 ) + (( edge_center[1] - (edge_length / 2.0) ) / edge_length) + 0.2 );

                    // unit outer normal creates the correct sign!
                    for( int k = 0; k < dimension; ++k )
                       jump_up_down[ id ] += pow( edge_length , 2.0 ) * ( diffusive_flux_y_S[ 0 ][ k ] * unitOuterNormal[ k ] );

                    //std :: cout << "edge_center = " << edge_center << std :: endl;
                    //std :: cout << "edge_length = " << edge_length << std :: endl;
                    //std :: cout << "unitOuterNormal = " << unitOuterNormal << std :: endl;
                    //std :: cout << "diffusive_flux_y_S = " << diffusive_flux_y_S[ 0 ] << std :: endl;
                    //std :: cout << "jump_up_down id = " << id << std :: endl;
                    //std :: cout << "jump_up_down[" << id << "] = " << jump_up_down[ id ] << std :: endl << std :: endl;

                    eps_delta_boundary_edge_found = true;
                    //std :: cout << "Found eps/delta boundary edge with center = " << edge_center << std :: endl;

                   }

                  if ( (fabs(edge_center[1]) == ((0.5*(epsilon_estimated/delta)))) &&
                       (fabs(edge_center[0]) < ((0.5*(epsilon_estimated/delta))))   )
                   {

                    // + 0.2 to avoid rounding errors!
                    int id = int( ( num_boundary_faces_per_direction / 2 ) + (( edge_center[0] - (edge_length / 2.0) ) / edge_length) + 0.2 );

                    // unit outer normal creates the correct sign!
                    for( int k = 0; k < dimension; ++k )
                       jump_left_right[ id ] += pow( edge_length , 2.0 ) * ( diffusive_flux_y_S[ 0 ][ k ] * unitOuterNormal[ k ] );

                    //std :: cout << "edge_center = " << edge_center << std :: endl;
                    //std :: cout << "edge_length = " << edge_length << std :: endl;
                    //std :: cout << "unitOuterNormal = " << unitOuterNormal << std :: endl;
                    //std :: cout << "diffusive_flux_y_S = " << diffusive_flux_y_S[ 0 ] << std :: endl;
                    //std :: cout << "jump_left_right id = " << id << std :: endl;
                    //std :: cout << "jump_left_right[" << id << "] = " << jump_left_right[ id ] << std :: endl << std :: endl;

                    eps_delta_boundary_edge_found = true;
                    //std :: cout << "Found eps/delta boundary edge with center = " << edge_center << std :: endl;

                   }

                }

               if ( ( fabs(edge_center[0]) < (0.5*(epsilon_estimated/delta)) ) && 
                    ( fabs(edge_center[1]) < (0.5*(epsilon_estimated/delta)) )  )
                {
                   // in this situation, p_nit always has an outside neighbor entity
                   // (use the normal gradient jumps)


                   // --------- the 'outer entity' (micro grid) ------------
                   //                ( neighbor entity )

                   EntityPointerType outer_p_it = p_nit->outside();
                   const EntityType &outer_micro_entity = *outer_p_it;

                   // one point quadrature formula ( A_h^{\eps}(y)=A^{\eps}(y_s) )
                   EntityQuadratureType outer_one_point_quadrature( outer_micro_entity , 0 );

                   // Q_h(u_H)(x_T,y) on the neighbor entity:
                   PeriodicLocalFunctionType outer_loc_Q_u_H_x_T = corrector_u_H_on_entity.localFunction(outer_micro_entity);
                   JacobianRangeType gradient_outer_Q_u_H_x_T(0.);
                   outer_loc_Q_u_H_x_T.jacobian(outer_one_point_quadrature[ 0 ], gradient_outer_Q_u_H_x_T);

                   // S denotes the micro grid element (i.e. 'micro_entity')
                   const EntityGeometryType& outer_geometry_S = outer_micro_entity.geometry();
                   // outer_y_S denotes the barycenter of the neighbor micro grid element of S:
                   const DomainType &outer_y_S = outer_geometry_S.global(outer_one_point_quadrature.point(0));

                   // to evaluate A^{\epsilon}_h (in center of current outer entity):
                   DomainType outer_globalPoint;
                   JacobianRangeType direction_of_diffusion_outside;
                   for( int k = 0; k < dimension; ++k )
                     {
                      outer_globalPoint[ k ] = x_T[ k ] + (delta * outer_y_S[ k ]);
                      direction_of_diffusion_outside[ 0 ][ k ] = gradient_u_H[ 0 ][ k ] +  gradient_outer_Q_u_H_x_T[ 0 ][ k ];
                     }

                   JacobianRangeType diffusive_flux_y_S_outside;
                   diffusion_.diffusiveFlux( outer_globalPoint,
                                             direction_of_diffusion_outside,
                                             diffusive_flux_y_S_outside );

                   // ----------------------------------------------------

                   RangeType aux_value = 0.0;

                   for( int k = 0; k < dimension; ++k )
                     aux_value += ( diffusive_flux_y_S_outside[ 0 ][ k ] - diffusive_flux_y_S[ 0 ][ k ] ) * unitOuterNormal[ k ] ;

                   local_indicator += pow( edge_length , 4.0 ) * pow( aux_value , 2.0 );
                }

            } // end Intersection iterator
        } // endif: "y_S \in \eps/\delta Y"

        }

       for( int id = 0; id < num_boundary_faces_per_direction; ++id )
         {
          local_indicator += ( pow( jump_up_down[ id ] , 2.0 ) + pow( jump_left_right[ id ] , 2.0 ) );
         }

       if ( eps_delta_boundary_edge_found == false )
        {
           std :: cout << "Error! Make sure that the restriction of the Y-triangulation on 'eps/delta Y' is a complete periodic triangulation of 'eps/delta Y' on its own. (for instance: delta = 2 epsilon should work)" << std :: endl; 
           std :: abort();
        }

       local_indicator *= entityVolume;

       return local_indicator;

     } // end of method
#endif

#endif

    // adaptive_refinement
    void adaptive_refinement( GridType &coarse_grid,
                              const DiscreteFunctionType& msfem_solution,
                              const DiscreteFunctionType& msfem_coarse_part,
                              const DiscreteFunctionType& msfem_fine_part,
                              std :: ofstream& data_file )
    {
      
       std :: cout << "Start computing conservative fluxes..." << std :: endl;
      
       ConservativeFluxProblemSolver< SubGridDiscreteFunctionType, DiscreteFunctionType, DiffusionOperatorType, MacroMicroGridSpecifierImp >
           flux_problem_solver( fineDiscreteFunctionSpace_, diffusion_, specifier_, data_file, path_ );

       flux_problem_solver.solve_all( subgrid_list_ );
       
       std :: cout << "Conservative fluxes computed successfully." << std :: endl;
       
      
       std :: cout << "Starting error estimation..." << std :: endl;

       const int number_of_coarse_grid_entities = specifier_.getNumOfCoarseEntities();
       const int number_of_fine_grid_entities = msfem_solution.space().gridPart().grid().size( 0 /*codim*/ );

       // error contribution of H^2 ||f||_L2(T):
       std :: vector < RangeType > loc_coarse_residual( number_of_coarse_grid_entities );

       std :: vector < RangeType > loc_projection_error( number_of_coarse_grid_entities );

       std :: vector < RangeType > loc_coarse_grid_jumps( number_of_coarse_grid_entities );

       std :: vector < RangeType > loc_conservative_flux_jumps( number_of_coarse_grid_entities );

       // Summation ueber die einzelnen fine-grid indicators die zu coarse grid entity beitragen:

       std :: vector < RangeType > loc_approximation_error( number_of_coarse_grid_entities );

       std :: vector < RangeType > loc_fine_grid_jumps( number_of_coarse_grid_entities );

       RangeType total_coarse_residual( 0.0 );
       RangeType total_projection_error( 0.0 );
       RangeType total_coarse_grid_jumps( 0.0 );
       RangeType total_conservative_flux_jumps( 0.0 );
       RangeType total_approximation_error( 0.0 );
       RangeType total_fine_grid_jumps( 0.0 );

       RangeType total_estimated_error( 0.0 );

       const DiscreteFunctionSpaceType& coarseDiscreteFunctionSpace = specifier_.coarseSpace();
       const LeafIndexSetType& coarseGridLeafIndexSet = coarseDiscreteFunctionSpace.gridPart().grid().leafIndexSet();

       // Coarse Entity Iterator 
       const IteratorType coarse_grid_end = coarseDiscreteFunctionSpace.end();
       for( IteratorType coarse_grid_it = coarseDiscreteFunctionSpace.begin(); coarse_grid_it != coarse_grid_end; ++coarse_grid_it )
        {
          int global_index_entity = coarseGridLeafIndexSet.index( *coarse_grid_it );

          loc_coarse_residual[ global_index_entity ] = indicator_f( *coarse_grid_it );
          total_coarse_residual += pow( loc_coarse_residual[ global_index_entity ], 2.0 );

	  getFluxes( *coarse_grid_it,
                     msfem_coarse_part,
                     loc_conservative_flux_jumps[ global_index_entity ],
                     loc_coarse_grid_jumps[ global_index_entity ] );

          total_conservative_flux_jumps += pow( loc_conservative_flux_jumps[ global_index_entity ], 2.0 );
          total_coarse_grid_jumps += pow( loc_coarse_grid_jumps[ global_index_entity ], 2.0 );

#if 1

          // the sub grid U(T) that belongs to the coarse_grid_entity T
          SubGridType& sub_grid_U_T = subgrid_list_.getSubGrid( global_index_entity );
          SubGridPart subGridPart( sub_grid_U_T );

          SubGridDiscreteFunctionSpaceType localDiscreteFunctionSpace( subGridPart );

          SubGridDiscreteFunctionType local_problem_solution_e0( "Local problem Solution e_0", localDiscreteFunctionSpace );
          local_problem_solution_e0.clear();

          SubGridDiscreteFunctionType local_problem_solution_e1( "Local problem Solution e_1", localDiscreteFunctionSpace );
          local_problem_solution_e1.clear();

          // --------- load local solutions -------

          char location_lps[50];
          sprintf( location_lps, "/local_problems/_localProblemSolutions_%d", global_index_entity );
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
          else
           { std :: cout << "Error! Could not read data file for the local problem solutions." << std :: endl; abort(); }

          if (reader_is_open)
           { discrete_function_reader.read( 1, local_problem_solution_e1 ); }

#endif

#if 1

          // iterator for the local micro grid ('the subgrid corresponding with U(T)')
          const SubGridIteratorType local_grid_it_end = localDiscreteFunctionSpace.end();
          for( SubGridIteratorType local_grid_it = localDiscreteFunctionSpace.begin(); local_grid_it != local_grid_it_end; ++local_grid_it )
              {

                const SubGridEntityType &local_grid_entity = *local_grid_it;

                // check if "local_grid_entity" (which is an entity of U(T)) is in T:
                // -------------------------------------------------------------------

                const EntityPointerType host_local_grid_it = localDiscreteFunctionSpace.grid().template getHostEntity<0>( local_grid_entity );

                EntityPointerType father_of_loc_grid_it = host_local_grid_it;

                for (int lev = 0; lev < specifier_.getLevelDifference() ; ++lev)
                       father_of_loc_grid_it = father_of_loc_grid_it->father();

                bool father_found = coarseGridLeafIndexSet.contains( *father_of_loc_grid_it );
                while ( father_found == false )
                 {
                   father_of_loc_grid_it = father_of_loc_grid_it->father();
                   father_found = coarseGridLeafIndexSet.contains( *father_of_loc_grid_it );
                 }

                bool entities_identical = true;
                int number_of_nodes = (*coarse_grid_it).template count<2>();
                for ( int k = 0; k < number_of_nodes; k += 1 )
                  {
                    if ( !(coarse_grid_it->geometry().corner(k) == father_of_loc_grid_it->geometry().corner(k)) )
                      { entities_identical = false; }
                  }

                if ( entities_identical == false )
                  continue;

                // -------------------------------------------------------------------

                const SubGridEntityGeometryType &local_grid_geometry = local_grid_entity.geometry();
                assert( local_grid_entity.partitionType() == InteriorEntity );
#if 0
		

		












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
#endif
              }
#endif
        }




// fine-grid iterator:
#if 1


       // Coarse Entity Iterator 
       const IteratorType fine_grid_end = fineDiscreteFunctionSpace_.end();
       for( IteratorType fine_grid_it = fineDiscreteFunctionSpace_.begin(); fine_grid_it != fine_grid_end; ++fine_grid_it )
        {

          EntityType& entity = *fine_grid_it;

          // identify coarse grid father entity
          EntityPointerType coarse_father = fine_grid_it;
          for (int lev = 0; lev < specifier_.getLevelDifference() ; ++lev)
            coarse_father = coarse_father->father();

          bool father_found = coarseGridLeafIndexSet.contains( *coarse_father );
          while ( father_found == false )
               {
                 coarse_father = coarse_father->father();
                 father_found = coarseGridLeafIndexSet.contains( *coarse_father );
               }

          int coarse_father_index = coarseGridLeafIndexSet.index( *coarse_father );

          const EntityGeometryType& entityGeometry = entity.geometry();

          EntityQuadratureType entityQuadrature( entity , 0 ); // 0 = polynomial order
          const DomainType &x = entityGeometry.global( entityQuadrature.point(0) );

          LocalFunctionType local_msfem_sol = msfem_solution.localFunction( entity );
          JacobianRangeType gradient_msfem_sol(0.);
          local_msfem_sol.jacobian( entityQuadrature[ 0 ], gradient_msfem_sol);

          JacobianRangeType diffusive_flux_x;
          diffusion_.diffusiveFlux( x, gradient_msfem_sol, diffusive_flux_x );

          EntityQuadratureType highOrder_entityQuadrature( entity , 2*spacePolOrd+2 );

          const int quadratureNop = highOrder_entityQuadrature.nop();
          for( int quadraturePoint = 0; quadraturePoint < quadratureNop; ++quadraturePoint )
            {
              const double weight = highOrder_entityQuadrature.weight( quadraturePoint ) *
                entityGeometry.integrationElement( highOrder_entityQuadrature.point( quadraturePoint ) );

              DomainType point = entityGeometry.global( highOrder_entityQuadrature.point( quadraturePoint ) );

              JacobianRangeType diffusive_flux_high_order;
              diffusion_.diffusiveFlux( point, gradient_msfem_sol, diffusive_flux_high_order );

              RangeType value = 0.0;
              for ( int i = 0; i < dimension; ++i )
                value += pow( diffusive_flux_x[ 0 ][ i ] - diffusive_flux_high_order[ 0 ][ i ], 2.0 );

              loc_approximation_error[ coarse_father_index ] += weight * value;
              total_approximation_error += weight * value;

            }


          const GridPartType &fineGridPart = fineDiscreteFunctionSpace_.gridPart();

          IntersectionIteratorType endnit = fineGridPart.iend( *fine_grid_it );
          for( IntersectionIteratorType nit = fineGridPart.ibegin( *fine_grid_it ); nit != endnit ; ++nit)
            {
              FaceQuadratureType innerFaceQuadrature( fineGridPart, *nit, 0 , FaceQuadratureType::INSIDE);
              const FaceGeometryType& faceGeometry = nit->geometry();

              if ( nit->neighbor() == false )
               { continue; }

              EntityPointerType outer_fine_grid_it = nit->outside();
              EntityType& outer_entity = *outer_fine_grid_it;

              EntityQuadratureType outer_entityQuadrature( outer_entity , 0 ); // 0 = polynomial order
              const EntityGeometryType& outer_entityGeometry = outer_entity.geometry();
              const DomainType &outer_x = outer_entityGeometry.global( outer_entityQuadrature.point(0) );

              LocalFunctionType outer_local_msfem_sol = msfem_solution.localFunction( outer_entity );
              JacobianRangeType outer_gradient_msfem_sol(0.);
              outer_local_msfem_sol.jacobian( outer_entityQuadrature[ 0 ], outer_gradient_msfem_sol);

              JacobianRangeType diffusive_flux_outside;
              diffusion_.diffusiveFlux( outer_x, outer_gradient_msfem_sol, diffusive_flux_outside );

              DomainType unitOuterNormal =
                 nit->unitOuterNormal( innerFaceQuadrature.localPoint( 0 ) );

              const RangeType edge_length = faceGeometry.volume();

              RangeType int_value = 0.0;
              for ( int i = 0; i < dimension; ++i )
                int_value[ i ] += (diffusive_flux_x[ 0 ][ i ] - diffusive_flux_outside[ 0 ][ i ]) * unitOuterNormal[ i ];
              int_value = pow( int_value, 2.0 );

              loc_fine_grid_jumps[ coarse_father_index ] += edge_length * edge_length * int_value;
              total_fine_grid_jumps += edge_length * edge_length * int_value;
            }

        }

    for ( int m = 0; m < number_of_coarse_grid_entities; ++m )
      {
        loc_approximation_error[ m ] = sqrt( loc_approximation_error[ m ] );
        loc_fine_grid_jumps[ m ] = sqrt( loc_fine_grid_jumps[ m ] );
      }

#endif








       total_coarse_residual = sqrt( total_coarse_residual );
       total_projection_error = sqrt( total_projection_error );
       total_coarse_grid_jumps = sqrt( total_coarse_grid_jumps );
       total_conservative_flux_jumps = sqrt( total_conservative_flux_jumps );
       total_approximation_error = sqrt( total_approximation_error );
       total_fine_grid_jumps = sqrt( total_fine_grid_jumps );

       total_estimated_error += total_coarse_residual;
       total_estimated_error += total_projection_error;
       total_estimated_error += total_coarse_grid_jumps;
       total_estimated_error += total_conservative_flux_jumps;
       total_estimated_error += total_approximation_error;
       total_estimated_error += total_fine_grid_jumps;

       if (data_file.is_open())
        {
           data_file << std :: endl;
           data_file << "Estimated Errors:" << std :: endl << std :: endl;
           data_file << "Total estimated error = " << total_estimated_error << "." << std :: endl;
           data_file << "where: " << std :: endl;
           data_file << "total_coarse_residual = " << total_coarse_residual << "." << std :: endl;
           data_file << "total_projection_error = " << total_projection_error << "." << std :: endl;
           data_file << "total_coarse_grid_jumps = " << total_coarse_grid_jumps << "." << std :: endl;
           data_file << "total_conservative_flux_jumps = " << total_conservative_flux_jumps << "." << std :: endl;
           data_file << "total_approximation_error = " << total_approximation_error << "." << std :: endl;
           data_file << "total_fine_grid_jumps = " << total_fine_grid_jumps << "." << std :: endl;
        }

       std :: cout << std :: endl;
       std :: cout << "Estimated Errors:" << std :: endl << std :: endl;
       std :: cout << "Total estimated error = " << total_estimated_error << "." << std :: endl;
       std :: cout << "where: " << std :: endl;
       std :: cout << "total_coarse_residual = " << total_coarse_residual << "." << std :: endl;
       std :: cout << "total_projection_error = " << total_projection_error << "." << std :: endl;
       std :: cout << "total_coarse_grid_jumps = " << total_coarse_grid_jumps << "." << std :: endl;
       std :: cout << "total_conservative_flux_jumps = " << total_conservative_flux_jumps << "." << std :: endl;
       std :: cout << "total_approximation_error = " << total_approximation_error << "." << std :: endl;
       std :: cout << "total_fine_grid_jumps = " << total_fine_grid_jumps << "." << std :: endl;

    }

}; // end of class MsFEMErrorEstimator

} // end namespace 

#endif
