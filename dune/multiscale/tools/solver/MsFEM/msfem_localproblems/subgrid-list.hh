#ifndef SUBGRIDLIST_HH
#define SUBGRIDLIST_HH

#include <dune/common/fmatrix.hh>

#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/multiscale/tools/solver/MsFEM/msfem_localproblems/localproblemsolver.hh>

#include <dune/fem/operator/2order/lagrangematrixsetup.hh>

namespace Dune
{

  template< class HostDiscreteFunctionImp, class SubGridImp >
  class SubGridList
  {
  public:
    //! ---------------- typedefs for the HostDiscreteFunctionSpace -----------------------

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
    
    typedef typename HostDiscreteFunctionSpaceType :: IteratorType MaxLevelHostIteratorType;
    
    typedef typename MaxLevelHostIteratorType :: Entity HostEntityType;

    typedef typename HostEntityType :: EntityPointer HostEntityPointerType;
    
    typedef typename HostGridType :: template Codim< 0 > :: template Partition< All_Partition > :: LevelIterator HostgridLevelEntityIteratorType;
   
    typedef typename HostGridType :: Traits :: LevelIndexSet HostGridLevelIndexSet;

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
		     const int& level_difference,
		     int &father_index,
                     const HostGridPartType& hostGridPart,
                     SubGridType& subGrid,
                     EntityPointerCollectionType& entities_sharing_same_node,
                     int& layer,
		     bool***& enriched )
    {

      for ( int l = 0; l <= layer; l+=1 )
         enriched[ father_index ][ hostGridPart.indexSet().index( *hit ) ][ l ] = true; 
   
      layer -= 1;
      // loop over the nodes of the enity
      for ( int i = 0; i < (*hit).template count<2>(); i += 1 )
	{

	   const HostNodePointer node = (*hit).template subEntity<2>(i);
	   int global_index_node = hostGridPart.indexSet().index( *node );

	   for( int j = 0; j < entities_sharing_same_node[global_index_node].size(); j += 1 )
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

		      if ( !( father == level_father_it ) )
		       { 
			 if ( enriched[ father_index ][ hostGridPart.indexSet().index( *entities_sharing_same_node[ global_index_node ][ j ] ) ][ layer ] == false )
			  {
			    enrichment( entities_sharing_same_node[ global_index_node ][ j ], level_father_it,
			    	        level_difference, father_index,
			                hostGridPart, subGrid, entities_sharing_same_node, layer, enriched );

     	                    layer += 1;
			  }
			 
		      }
		   }

	      }

	}

   }
    
    SubGridList( const HostDiscreteFunctionSpaceType& hostSpace, const std :: vector < int >& number_of_layers, const int& computational_level, bool silent = true )
    : hostSpace_( hostSpace ),
      number_of_layers_( number_of_layers ),
      computational_level_( computational_level ),
      silent_( silent )
    {
      
      const HostGridPartType& hostGridPart = hostSpace_.gridPart();

      HostGridType& hostGrid = hostSpace_.gridPart().grid();

      int number_of_nodes = hostGrid.size( hostGrid.maxLevel(), 2 /*codim*/ );
      
      // -------- identify the entities that share a certain node -------

      std :: vector< std :: vector < HostEntityPointerType > > entities_sharing_same_node( number_of_nodes );
      
      for( MaxLevelHostIteratorType it = hostSpace_.begin(); it != hostSpace_.end(); ++it )
        {
	  int number_of_nodes_in_entity = (*it).template count<2>();
	  for ( int i = 0; i < number_of_nodes_in_entity; i += 1 )
	    {
	      const HostNodePointer node = (*it).template subEntity<2>(i);
	      int global_index_node = hostGridPart.indexSet().index( *node );
	      
	      entities_sharing_same_node[ global_index_node ].push_back( it );
	    }
        }

      int max_num_layers = 0;
      for ( int i = 0; i < number_of_layers.size(); i += 1 )
        { if ( number_of_layers[i] > max_num_layers )
	   { max_num_layers = number_of_layers[i]; }
	}
        
#if 0
      for ( int i = 0; i < number_of_nodes ; i+=1 )
        {
          std :: cout << "Knoten " << i << " wird von " << entities_sharing_same_node[i].size() << " Entities geteilt." << std :: endl;
        }
#endif

      // ---------------------------------------------------------------- 
     
      
      // maximum level defined in this grid. Levels are numbered 0 ... maxLevel with 0 the coarsest level.
      int maxLevel = hostGrid.maxLevel();
      int level_difference = maxLevel - computational_level_;

      // number of grid entities of a given codim on a given level in this process.
      int number_of_level_host_entities = hostGrid.size( computational_level_, 0 /*codim*/ );


      std :: cout << "number_of_level_host_entities = " << number_of_level_host_entities << std :: endl;
      

      subGrid = new SubGridType* [ number_of_level_host_entities ];
      
      //! ----------- create subgrids --------------------

      const HostGridLevelIndexSet& hostGridLevelIndexSet = hostGrid.levelIndexSet(computational_level_);

      //SubGridType* subGrid[ number_of_level_host_entities ];

      HostgridLevelEntityIteratorType level_iterator_end = hostGrid.template lend< 0 >( computational_level_ );
      HostgridLevelEntityIteratorType level_iterator_begin = hostGrid.template lbegin< 0 >( computational_level_ );

      int lev_index = 0;
      for ( HostgridLevelEntityIteratorType lit = level_iterator_begin;
         lit != level_iterator_end ; ++lit )
        {

          subGrid[ lev_index ] = new SubGridType( hostGrid );
          subGrid[ lev_index ]->createBegin();

          //level_entity_collection.push_back( lit );
          lev_index += 1;

        }

      for( MaxLevelHostIteratorType it = hostSpace_.begin(); it != hostSpace_.end(); ++it )
        {
	  int number_of_nodes_in_entity = (*it).template count<2>();
	  for ( int i = 0; i < number_of_nodes_in_entity; i += 1 )
	    {
	      const HostNodePointer node = (*it).template subEntity<2>(i);
	      int global_index_node = hostGridPart.indexSet().index( *node );
	      
	      entities_sharing_same_node[ global_index_node ].push_back( it );
	    }
        }
      
       bool *** enriched = new bool ** [ number_of_level_host_entities ];
       for ( int k = 0; k < number_of_level_host_entities ; k += 1 )
       {
	 enriched[ k ] = new bool * [ hostGrid.size( hostGrid.maxLevel(), 0 ) ];
         for ( int m = 0; m < hostGrid.size( hostGrid.maxLevel(), 0 ) ; m += 1 )
	 {
	    enriched[ k ][ m ] = new bool [ max_num_layers + 1 ];
	 }
       }
       
       for ( int k = 0; k < number_of_level_host_entities ; k += 1 )
       {
         for ( int m = 0; m < hostGrid.size( hostGrid.maxLevel(), 0 ) ; m += 1 )
	 {
	    for ( int l = 0; l < max_num_layers + 1 ; l += 1 )
	    enriched[ k ][ m ][ l ] = false;
	 }
       }

	   
	   
      // a maxlevel iterator for the codim 0 hostgrid entities:
      MaxLevelHostIteratorType host_endit = hostSpace_.end();
      for( MaxLevelHostIteratorType host_it = hostSpace_.begin();
           host_it != host_endit ;
           ++host_it )
        {

           const HostEntityType& host_entity = *host_it;
	   
           int number_of_nodes_in_entity = (*host_it).template count<2>();
	   

           // get the level 'computational_level'-father of host_entity (which is a maxlevel entity)
           HostEntityPointerType level_father_entity = host_it;
           for (int lev = 0; lev < level_difference; ++lev)
             level_father_entity = level_father_entity->father();

           int father_index = hostGridLevelIndexSet.index( *level_father_entity );
           // std :: cout << "father_index = " << father_index << std :: endl;
	   

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
                   const HostEntityType& neighborHostEntity = *neighborHostEntityPointer;
		   
                   HostEntityPointerType level_father_neighbor_entity = neighborHostEntityPointer;
                   for (int lev = 0; lev < level_difference; ++lev)
                       level_father_neighbor_entity = level_father_neighbor_entity->father();
	   
	   
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

	   int layers = number_of_layers[ father_index ];
	   if ( layers > 0 )
	    {
	       enrichment( host_it, level_father_entity, level_difference, father_index,
		           hostGridPart, *subGrid[ father_index ], entities_sharing_same_node, layers, enriched );
	    }
            
        }

      for ( int i = 0; i < number_of_level_host_entities; ++i )
       {
          subGrid[ i ]->createEnd();
          if ( !silent_ )
            { subGrid[ i ]->report(); }
       }

      //! ----------- end create subgrids --------------------

    }
    
  SubGridType& getSubGrid( int i )
  {
    int size = hostSpace_.gridPart().grid().size( computational_level_, 0 /*codim*/ );
    
    if ( i < size )    
     { return *(subGrid[ i ]); }
    else
     { std :: cout << "Error. Subgrid-Index too large." << std :: endl; }
  }
    
  private:
    
    const HostDiscreteFunctionSpaceType &hostSpace_;
    const std :: vector < int >& number_of_layers_;
    const int computational_level_;
    
    bool silent_;
    
    SubGridType** subGrid;
    
  public:
    
    //Destruktor (Dynamisch angeforderter Speicher wird wieder freigegeben)
    ~SubGridList()
     {
        int size = hostSpace_.gridPart().grid().size( computational_level_, 0 /*codim*/ );
        
        for (unsigned int i = 0;i<size; ++i)
         delete subGrid[i];
        delete[] subGrid;

     }

  };

}

#endif // #ifndef SUBGRIDLIST_HH
