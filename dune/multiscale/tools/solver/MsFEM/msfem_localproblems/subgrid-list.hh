#ifndef SUBGRIDLIST_HH
#define SUBGRIDLIST_HH

#include <dune/common/fmatrix.hh>
#include <dune/subgrid/subgrid.hh>

#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/multiscale/tools/solver/MsFEM/msfem_localproblems/localproblemsolver.hh>

#include <dune/fem/operator/2order/lagrangematrixsetup.hh>

/// done

namespace Dune
{

  template< class HostDiscreteFunctionImp, class SubGridImp, class MacroMicroGridSpecifierImp >
  class SubGridList
  {
  public:
    //! ---------------- typedefs for the HostDiscreteFunctionSpace -----------------------

    typedef MacroMicroGridSpecifierImp MacroMicroGridSpecifierType;
    
    typedef HostDiscreteFunctionImp HostDiscreteFunctionType;
    
    //! type of discrete function space
    typedef typename HostDiscreteFunctionType :: DiscreteFunctionSpaceType
      HostDiscreteFunctionSpaceType;

    //! type of (non-discrete )function space
    typedef typename HostDiscreteFunctionSpaceType :: FunctionSpaceType FunctionSpaceType;

    //! type of grid partition
    typedef typename HostDiscreteFunctionSpaceType :: GridPartType HostGridPartType;

    //! type of grid
    typedef typename HostDiscreteFunctionSpaceType :: GridType HostGridType;

    typedef typename HostGridType :: Traits :: LeafIndexSet HostGridLeafIndexSet;

    typedef typename HostDiscreteFunctionSpaceType :: IteratorType HostGridEntityIteratorType;

    typedef typename HostGridEntityIteratorType :: Entity HostEntityType;

    typedef typename HostEntityType :: EntityPointer HostEntityPointerType;

// old:
#if 0
    typedef typename HostGridType :: template Codim< 0 > :: template Partition< All_Partition > :: LevelIterator HostgridLevelEntityIteratorType;
    typedef typename HostGridType :: Traits :: LevelIndexSet HostGridLevelIndexSet;

    //usage:
    // const HostGridLevelIndexSet& hostGridLevelIndexSet = hostGrid.levelIndexSet(computational_level_);
    // int father_index = hostGridLevelIndexSet.index( *level_it );
#endif

    typedef typename HostEntityType :: template Codim< 2 > :: EntityPointer HostNodePointer;

    typedef typename HostGridPartType :: IntersectionIteratorType HostIntersectionIterator;

    //! ---------------- typedefs for the SubgridDiscreteFunctionSpace -----------------------
    //  ( typedefs for the local grid and the corresponding local ('sub') )discrete space ) 

    //! type of grid
    typedef SubGridImp SubGridType; 

    //! type of grid part
    typedef LeafGridPart< SubGridType > SubGridPartType;

    
    template< typename EntityPointerCollectionType >
    bool entity_patch_in_subgrid( HostEntityPointerType& hit,
                           const HostGridPartType& hostGridPart,
                           SubGridType& subGrid,
                           EntityPointerCollectionType& entities_sharing_same_node )
    {

      bool patch_in_subgrid_ = true;
      // loop over the nodes of the enity
      for ( int i = 0; i < (*hit).template count<2>(); i += 1 )
	{

	   const HostNodePointer node = (*hit).template subEntity<2>(i);

	   int global_index_node = hostGridPart.indexSet().index( *node );

	   for( int j = 0; j < entities_sharing_same_node[global_index_node].size(); j += 1 )
	      {

		 if ( !( subGrid.template contains <0>( *entities_sharing_same_node[ global_index_node ][ j ] ) ) )
		   {
		      patch_in_subgrid_ = false;
		   }

	      }

	}
      return patch_in_subgrid_;

   }    



    template< typename EntityPointerCollectionType >
    void enrichment( HostEntityPointerType& hit,
		     HostEntityPointerType& level_father_it,
		     MacroMicroGridSpecifierType& specifier,
		     int &father_index,
                     const HostGridPartType& hostGridPart,
                     SubGridType& subGrid,
                     EntityPointerCollectionType& entities_sharing_same_node,
                     int& layer,
		     bool***& enriched )
    {

      // difference in levels between coarse and fine grid
      int level_difference = specifier.getLevelDifference();
      HostDiscreteFunctionSpaceType& coarseSpace = specifier.coarseSpace();

      const HostGridLeafIndexSet& coarseGridLeafIndexSet = coarseSpace.gridPart().grid().leafIndexSet();

      const HostGridLeafIndexSet& hostGridLeafIndexSet = hostSpace_.gridPart().grid().leafIndexSet();

      for ( int l = 0; l <= layer; l+=1 )
        {

          enriched[ father_index ][ hostGridLeafIndexSet.index( *hit ) ][ l ] = true;

         }

      layer -= 1;


      // loop over the nodes of the enity
      for ( int i = 0; i < (*hit).template count<2>(); i += 1 )
	{

           const HostNodePointer node = (*hit).template subEntity<2>(i);
	   int global_index_node = hostGridPart.indexSet().index( *node );

       for( size_t j = 0; j < entities_sharing_same_node[global_index_node].size(); j += 1 )
	     {

               if ( !( subGrid.template contains <0>( *entities_sharing_same_node[ global_index_node ][ j ] ) ) )
		{
		  subGrid.insertPartial( *entities_sharing_same_node[ global_index_node ][ j ] );
                }

               if ( layer > 0 )
                {

                  HostEntityPointerType father = entities_sharing_same_node[ global_index_node ][ j ];
                  for (int lev = 0; lev < level_difference; ++lev)
                       father = father->father();

//! new version:
#if 1
                  HostEntityPointerType coarse_father_test = father;
	
                  bool father_found = false;
                  while ( father_found == false )
                    {

                      if ( coarseGridLeafIndexSet.contains( *coarse_father_test ) == true )
                       {
                         father = coarse_father_test;
                       }

                     if ( coarse_father_test->hasFather() == false )
                       { father_found = true; }
                     else
                       { coarse_father_test = coarse_father_test->father(); }
	        
                    }
#endif

//! old version
#if 0
                  bool father_found = false;
                  while ( father_found == false )
                   {

                     if ( coarseGridLeafIndexSet.contains( *father ) == true )
                      {
                        if ( father->hasFather() == false )
                         { father_found = true; }
                        else
                         {
                           if ( coarseGridLeafIndexSet.contains( *(father->father())) == false )
                            { father_found = true; }
                         }
                      }

                     if ( father_found == false )
                        { father = father->father(); }

                   }
#endif

                  if ( !( father == level_father_it ) )
                   {

                     if ( enriched[ father_index ][ hostGridLeafIndexSet.index( *entities_sharing_same_node[ global_index_node ][ j ] ) ][ layer ] == false )
                      {
                        enrichment( entities_sharing_same_node[ global_index_node ][ j ], level_father_it,
                        specifier, father_index,
                        hostGridPart, subGrid, entities_sharing_same_node, layer, enriched );

                        layer += 1;
                      }


                   }
                }
	      }
	}

   }
    
    SubGridList( MacroMicroGridSpecifierType& specifier, bool silent = true )
    : hostSpace_( specifier.fineSpace() ),
      specifier_( specifier ),
      silent_( silent )
    {

      std :: cout << "Starting creation of subgrids." << std :: endl << std :: endl; 

      HostDiscreteFunctionSpaceType& coarseSpace = specifier.coarseSpace();

      const HostGridPartType& hostGridPart = hostSpace_.gridPart();

      HostGridType& hostGrid = hostSpace_.gridPart().grid();

      int number_of_nodes = hostGrid.size( 2 /*codim*/ );

      // -------- identify the entities that share a certain node -------

      std :: vector< std :: vector < HostEntityPointerType > > entities_sharing_same_node( number_of_nodes );
      
      for( HostGridEntityIteratorType it = hostSpace_.begin(); it != hostSpace_.end(); ++it )
        {
	  int number_of_nodes_in_entity = (*it).template count<2>();
	  for ( int i = 0; i < number_of_nodes_in_entity; i += 1 )
	    {
	      const HostNodePointer node = (*it).template subEntity<2>(i);
	      int global_index_node = hostGridPart.indexSet().index( *node );

          entities_sharing_same_node[ global_index_node ].push_back( HostEntityPointerType(*it) );
	    }
        }

      int max_num_layers = 0;
      for ( int i = 0; i < specifier_.getNumOfCoarseEntities(); i += 1 )
        {
          if ( specifier_.getLayer(i) > max_num_layers )
           { max_num_layers = specifier_.getLayer(i); }
        }

      #if 0
      for ( int i = 0; i < number_of_nodes ; i+=1 )
        {
          std :: cout << "Node " << i << " is shared by " << entities_sharing_same_node[i].size() << " entities." << std :: endl;
        }
      #endif

      // ---------------------------------------------------------------- 
     
      
      // difference in levels between coarse and fine grid
      int level_difference = specifier.getLevelDifference();

      // number of coarse grid entities (of codim 0).
      int number_of_coarse_grid_entities = specifier.getNumOfCoarseEntities();


      std :: cout << "number_of_coarse_grid_entities = " << number_of_coarse_grid_entities << std :: endl;
      

      subGrid = new SubGridType* [ number_of_coarse_grid_entities ];
      
      //! ----------- create subgrids --------------------

      const HostGridLeafIndexSet& coarseGridLeafIndexSet = coarseSpace.gridPart().grid().leafIndexSet();

      for( HostGridEntityIteratorType coarse_it = coarseSpace.begin(); coarse_it != coarseSpace.end(); ++coarse_it )
        {

          int coarse_index = coarseGridLeafIndexSet.index( *coarse_it );

          subGrid[ coarse_index ] = new SubGridType( hostGrid );
          subGrid[ coarse_index ]->createBegin();
	}

      for( HostGridEntityIteratorType it = hostSpace_.begin(); it != hostSpace_.end(); ++it )
        {
	  int number_of_nodes_in_entity = (*it).template count<2>();
	  for ( int i = 0; i < number_of_nodes_in_entity; i += 1 )
	    {
	      const HostNodePointer node = (*it).template subEntity<2>(i);
	      int global_index_node = hostGridPart.indexSet().index( *node );

          entities_sharing_same_node[ global_index_node ].push_back( HostEntityPointerType(*it) );
	    }
        }
   
       bool *** enriched = new bool ** [ number_of_coarse_grid_entities ];
       for ( int k = 0; k < number_of_coarse_grid_entities ; k += 1 )
       {
	 enriched[ k ] = new bool * [ hostGrid.size( 0 /*codim*/ ) ];
         for ( int m = 0; m < hostGrid.size( 0 /*codim*/ ) ; m += 1 )
	 {
	    enriched[ k ][ m ] = new bool [ max_num_layers + 1 ];
	 }
       }

       for ( int k = 0; k < number_of_coarse_grid_entities ; k += 1 )
       {
         for ( int m = 0; m < hostGrid.size( 0 /*codim*/ ) ; m += 1 )
	 {
	    for ( int l = 0; l < max_num_layers + 1 ; l += 1 )
	    enriched[ k ][ m ][ l ] = false;
	 }
       }

      // a fine grid iterator for the codim 0 hostgrid entities:
      HostGridEntityIteratorType host_endit = hostSpace_.end();
      for( HostGridEntityIteratorType host_it = hostSpace_.begin();
           host_it != host_endit ;
           ++host_it )
        {

           const HostEntityType& host_entity = *host_it;

           int DUNE_UNUSED(number_of_nodes_in_entity) = (*host_it).template count<2>();

	   #if 0
           std :: cout << "host_it->geometry().corner(0) = " << host_it->geometry().corner(0) << std :: endl;
           std :: cout << "host_it->geometry().corner(1) = " << host_it->geometry().corner(1) << std :: endl;
           std :: cout << "host_it->geometry().corner(2) = " << host_it->geometry().corner(2) << std :: endl;
        #endif

           // get the coarse-grid-father of host_entity (which is a maxlevel entity)
           HostEntityPointerType level_father_entity = HostEntityPointerType(*host_it);



           for (int lev = 0; lev < level_difference; ++lev)
	     level_father_entity = level_father_entity->father();

           // changed 'contains'-method in 'indexset.hh'
           // we use: "return ( (subIndex >= 0) && (subIndex < size( codim )) );"
           // instead of "return (subIndex >= 0);"

	     
//! new version:
#if 1
           HostEntityPointerType coarse_father_test = level_father_entity;
	
           bool father_found = false;
           while ( father_found == false )
               {

                 if ( coarseGridLeafIndexSet.contains( *coarse_father_test ) == true )
                   { level_father_entity = coarse_father_test; }

                 if ( coarse_father_test->hasFather() == false )
                   { father_found = true; }
                 else
                   { coarse_father_test = coarse_father_test->father(); }
	        
               }
#endif

//! old version
#if 0
           bool father_found = false;
	   while ( father_found == false )
           {


             if ( coarseGridLeafIndexSet.contains( *level_father_entity ) == true )
              {
                  if ( level_father_entity->hasFather() == false )
                   { father_found = true; }
                  else
                   {
                     if ( coarseGridLeafIndexSet.contains( *(level_father_entity->father())) == false )
                       { father_found = true; }
                   }
               }

             if ( father_found == false )
               { level_father_entity = level_father_entity->father(); }
	   }
#endif

           #if 0
           std :: cout << "level_father_entity->geometry().corner(0) = " << level_father_entity->geometry().corner(0) << std :: endl;
           std :: cout << "level_father_entity->geometry().corner(1) = " << level_father_entity->geometry().corner(1) << std :: endl;
           std :: cout << "level_father_entity->geometry().corner(2) = " << level_father_entity->geometry().corner(2) << std :: endl << std :: endl;
           #endif

           int father_index = coarseGridLeafIndexSet.index( *level_father_entity );


           if ( !( subGrid[ father_index ]->template contains <0>(host_entity) ) )
            { subGrid[ father_index ]->insertPartial( host_entity ); }

           // check the neighbor entities and look if they belong to the same father
           // if yes, continue
           // if not, enrichement with 'n(T)'-layers
           bool all_neighbors_have_same_father = true;
           const HostIntersectionIterator iend = hostGridPart.iend( host_entity );
           for( HostIntersectionIterator iit = hostGridPart.ibegin( host_entity ); iit != iend; ++iit )
             {

                if ( iit->neighbor() ) //if there is a neighbor entity
                  {
                    // check if the neighbor entity is in the subgrid
                   const HostEntityPointerType neighborHostEntityPointer = iit->outside();
                   const HostEntityType& DUNE_UNUSED(neighborHostEntity) = *neighborHostEntityPointer;

                   HostEntityPointerType level_father_neighbor_entity = neighborHostEntityPointer;
                   for (int lev = 0; lev < level_difference; ++lev)
                       level_father_neighbor_entity = level_father_neighbor_entity->father();

//! new version:
#if 1
                   HostEntityPointerType coarse_father_test = level_father_neighbor_entity;
	
                   bool father_found = false;
                   while ( father_found == false )
                     {

                       if ( coarseGridLeafIndexSet.contains( *coarse_father_test ) == true )
                        { level_father_neighbor_entity = coarse_father_test; }

                       if ( coarse_father_test->hasFather() == false )
                        { father_found = true; }
                       else
                        { coarse_father_test = coarse_father_test->father(); }

                     }
#endif

//! old version
#if 0
	           father_found = false;
                   while ( father_found == false )
                    {

                      if ( coarseGridLeafIndexSet.contains( *level_father_neighbor_entity ) == true )
                       {
                         if ( level_father_neighbor_entity->hasFather() == false )
                          { father_found = true; }
                         else
                          {
                            if ( coarseGridLeafIndexSet.contains( *(level_father_neighbor_entity->father())) == false )
                              { father_found = true; }
                          }
                       }

                      if ( father_found == false )
                         { level_father_neighbor_entity = level_father_neighbor_entity->father(); }

                    }
#endif

                   if ( !(level_father_neighbor_entity == level_father_entity) )
                    {
                      all_neighbors_have_same_father = false;
                    }

                  }
                else
		  {
		    all_neighbors_have_same_father = false;
		  }
		
              }

           if ( all_neighbors_have_same_father == true )
	      { continue; }

	   int layers = specifier_.getLayer( father_index );
	   if ( layers > 0 )
	    {
           HostEntityPointerType hep(*host_it);
           enrichment( hep, level_father_entity, specifier, father_index,
		           hostGridPart, *subGrid[ father_index ], entities_sharing_same_node, layers, enriched );
	    }


        }

      for ( int i = 0; i < number_of_coarse_grid_entities; ++i )
       {
          subGrid[ i ]->createEnd();

#if 0
          if ( subGrid[ i ]->size( 2 ) == 0 )
            {
              subGrid[ i ]->createBegin();
              subGrid[ i ]->insertPartial( *hostSpace_.begin() );
              subGrid[ i ]->createEnd();
            }
#endif

          if ( !silent_ )
            {
              std :: cout << "Subgrid " << i << ":" << std :: endl;
              subGrid[ i ]->report();
            }

          if ( subGrid[ i ]->size( 2 ) == 0 )
            {
              std :: cout << "Error." << std :: endl;
              std :: cout << "Error. Created Subgrid with 0 nodes." << std :: endl;

              for( HostGridEntityIteratorType coarse_it = specifier_.coarseSpace().begin();
                   coarse_it != specifier_.coarseSpace().end() ; ++coarse_it )
                {

                   const HostEntityType& coarse_entity = *coarse_it;
                   int index = coarseGridLeafIndexSet.index( coarse_entity );
                   if ( i == index )
                    {
                      std :: cout << "We have a problem with the following coarse-grid element:" << std :: endl;
                      std :: cout << "coarse element corner(0) = " << coarse_it->geometry().corner(0) << std :: endl;
                      std :: cout << "coarse element corner(1) = " << coarse_it->geometry().corner(1) << std :: endl;
                      std :: cout << "coarse element corner(2) = " << coarse_it->geometry().corner(2) << std :: endl << std :: endl;
                    }
                }
              abort();
            }



       }

      //! ----------- end create subgrids --------------------

    }
    
// Kopierkonstruktor klappt nicht, da SubGrid keinen passenden Kopierkonstruktor besitzt
#if 0
  SubGridList( const SubGridList& list )
  : hostSpace_( list.hostSpace_ ),
    specifier_( list.specifier_ ),
    silent_( list.silent_ )
  {
      // number of coarse grid entities (of codim 0).
      int number_of_coarse_grid_entities = this->specifier_.getNumOfCoarseEntities();
      this->subGrid = new SubGridType* [ number_of_coarse_grid_entities ];
      
      for ( int i = 0; i < number_of_coarse_grid_entities; ++i )
       {
	  subGrid[ i ] = new SubGridType( *(list.subGrid[ i ]) );
       }
  }
#endif
    
  SubGridType& getSubGrid( int i )
  {
    int size = specifier_.getNumOfCoarseEntities();

    if ( i < size )
     { return *(subGrid[ i ]); }
    else
     { std :: cout << "Error. Subgrid-Index too large." << std :: endl; }
  }
    
  private:
    
    const HostDiscreteFunctionSpaceType &hostSpace_;
    MacroMicroGridSpecifierType& specifier_;
    
    bool silent_;
    
    SubGridType** subGrid;
    
  public:
    
    //Destruktor (Dynamisch angeforderter Speicher wird wieder freigegeben)
    ~SubGridList()
     {
        const int size = specifier_.getNumOfCoarseEntities();
        
        for ( int i = 0; i < size; ++i )
         delete subGrid[i];
        delete[] subGrid;

     }

  };

}

#endif // #ifndef SUBGRIDLIST_HH
