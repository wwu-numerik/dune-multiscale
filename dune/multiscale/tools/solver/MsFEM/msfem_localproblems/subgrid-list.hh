#ifndef SUBGRIDLIST_HH
#define SUBGRIDLIST_HH

#include <boost/noncopyable.hpp>
#include <boost/multi_array.hpp>

#include <dune/common/fmatrix.hh>
#include <dune/common/shared_ptr.hh>
#include <dune/common/exceptions.hh>
#include <dune/stuff/grid/entity.hh>
#include <dune/subgrid/subgrid.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/space/lagrangespace.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/operator/2order/lagrangematrixsetup.hh>
#include <dune/multiscale/tools/misc.hh>
#include <dune/multiscale/tools/subgrid_io.hh>
#include <dune/stuff/common/profiler.hh>
#include <dune/stuff/aliases.hh>
#include <dune/stuff/common/logging.hh>
#include <dune/stuff/common/parameter/configcontainer.hh>

namespace Dune {
//! container for cell problem subgrids
template< class HostDiscreteFunctionImp, class SubGridImp, class MacroMicroGridSpecifierImp >
class SubGridList : boost::noncopyable
{
public:
  //! ---------------- typedefs for the HostDiscreteFunctionSpace -----------------------

  typedef MacroMicroGridSpecifierImp MacroMicroGridSpecifierType;

  typedef HostDiscreteFunctionImp HostDiscreteFunctionType;

  //! type of discrete function space
  typedef typename HostDiscreteFunctionType::DiscreteFunctionSpaceType
  HostDiscreteFunctionSpaceType;

  //! type of (non-discrete )function space
  typedef typename HostDiscreteFunctionSpaceType::FunctionSpaceType FunctionSpaceType;

  //! type of domain
  typedef typename FunctionSpaceType::DomainType DomainType;
  
  //! type of grid partition
  typedef typename HostDiscreteFunctionSpaceType::GridPartType HostGridPartType;

  //! type of grid
  typedef typename HostDiscreteFunctionSpaceType::GridType HostGridType;

  typedef typename HostGridType::Traits::LeafIndexSet HostGridLeafIndexSet;

  typedef typename HostDiscreteFunctionSpaceType::IteratorType HostGridEntityIteratorType;

  typedef typename HostGridEntityIteratorType::Entity HostEntityType;

  typedef typename HostEntityType::EntityPointer HostEntityPointerType;

  typedef typename HostEntityType::template Codim< 2 >::EntityPointer HostNodePointer;

  typedef typename HostGridPartType::IntersectionIteratorType HostIntersectionIterator;

  typedef std::vector< DomainType> CoarseNodeVectorType;
  typedef std::vector< CoarseNodeVectorType > CoarseGridNodeStorageType;
  typedef boost::multi_array<bool, 3> EnrichmentMatrixType;
  
  //! ---------------- typedefs for the SubgridDiscreteFunctionSpace -----------------------
  // ( typedefs for the local grid and the corresponding local ('sub') )discrete space )

  //! type of grid
  typedef SubGridImp SubGridType;

  //! type of grid part
  typedef LeafGridPart< SubGridType > SubGridPartType;
  
    //! type of subgrid discrete function space
  typedef LagrangeDiscreteFunctionSpace< FunctionSpaceType, SubGridPartType, 1/*=POLORDER*/ > SubGridDiscreteFunctionSpace;

  //! type of subgrid discrete function
  typedef AdaptiveDiscreteFunction< SubGridDiscreteFunctionSpace > SubGridDiscreteFunction;
  
private:
  template< typename EntityPointerCollectionType >
  bool entity_patch_in_subgrid(const HostEntityPointerType& hit,
                               const HostGridPartType& hostGridPart,
                               shared_ptr<const SubGridType> subGrid,
                               const EntityPointerCollectionType& entities_sharing_same_node) const {
    bool patch_in_subgrid_ = true;

    // loop over the nodes of the enity
    for (int i = 0; i < (*hit).template count< 2 >(); i += 1)
    {
      const HostNodePointer node = (*hit).template subEntity< 2 >(i);

      const int global_index_node = hostGridPart.indexSet().index(*node);

      for (int j = 0; j < entities_sharing_same_node[global_index_node].size(); j += 1)
      {
        if ( !( subGrid.template contains< 0 >(*entities_sharing_same_node[global_index_node][j]) ) )
        {
          patch_in_subgrid_ = false;
        }
      }
    }
    return patch_in_subgrid_;
  } // entity_patch_in_subgrid

  /**
   * \note called in SubGridList constructor only
   */
  template< typename EntityPointerCollectionType >
  void enrichment(const HostEntityPointerType& hit,
                  const HostEntityPointerType& level_father_it,
                  const int& father_index, // father_index = index/number of current subgrid
                  const HostGridPartType& hostGridPart,
                  shared_ptr<SubGridType> subGrid,
                  EntityPointerCollectionType& entities_sharing_same_node,
                  int& layer,
                  EnrichmentMatrixType& enriched) {
    // difference in levels between coarse and fine grid
    const int level_difference = specifier_.getLevelDifference();
    HostDiscreteFunctionSpaceType& coarseSpace = specifier_.coarseSpace();

    const HostGridLeafIndexSet& coarseGridLeafIndexSet = coarseSpace.gridPart().grid().leafIndexSet();

    const HostGridLeafIndexSet& hostGridLeafIndexSet = hostSpace_.gridPart().grid().leafIndexSet();

    for (int l = 0; l <= layer; l += 1)
    {
      enriched[father_index][hostGridLeafIndexSet.index(*hit)][l] = true;
    }

    layer -= 1;

    // loop over the nodes of the fine grid enity
    for (int i = 0; i < (*hit).template count< 2 >(); i += 1)
    {
      const HostNodePointer node = (*hit).template subEntity< 2 >(i);
      int global_index_node = hostGridPart.indexSet().index(*node);

      for (size_t j = 0; j < entities_sharing_same_node[global_index_node].size(); j += 1)
      {
        if ( !( subGrid->template contains< 0 >(*entities_sharing_same_node[global_index_node][j]) ) )
        {
          subGrid->insertPartial(*entities_sharing_same_node[global_index_node][j]);

	  // get the corners of the father of the fine grid entity 'entities_sharing_same_node[global_index_node][j]'
	  // and add these corners to the vector 'coarse_node_store_[father_index]' (if they are not yet contained)
          if ( DSC_CONFIG_GET( "msfem.oversampling_strategy", 1 ) == 3 )
           {
              HostEntityPointerType coarse_father = Stuff::Grid::make_father(coarseGridLeafIndexSet,
                                                                             entities_sharing_same_node[global_index_node][j],
                                                                             level_difference);
              for (int c = 0; c < coarse_father->geometry().corners(); ++c )
	      {
                //! not an effective search algorithm (should be improved eventually):
                bool node_contained = false;
                for (size_t cn = 0; cn < coarse_node_store_[father_index].size(); ++cn )
	          {
                    // hard coding - 2d case:
                    if ( (coarse_node_store_[father_index][ cn ][ 0 ] == coarse_father->geometry().corner(c)[0])
		       && (coarse_node_store_[father_index][ cn ][ 1 ] == coarse_father->geometry().corner(c)[1] ) )
		      { node_contained = true; }
		  }
                if ( node_contained == false )
                 { coarse_node_store_[father_index].push_back( coarse_father->geometry().corner(c) ); }
	      }
            }

        }

        if (layer > 0)
        {
          HostEntityPointerType father = Stuff::Grid::make_father(coarseGridLeafIndexSet,
                                                                  entities_sharing_same_node[global_index_node][j],
                                                                  level_difference);
          if ( !(father == level_father_it) )
          {
            const auto& tmp_entity_ptr = entities_sharing_same_node[global_index_node][j];
            if (!enriched[father_index][hostGridLeafIndexSet
                .index(*tmp_entity_ptr)][layer])
            {
              enrichment(tmp_entity_ptr, level_father_it,
                         /*specifier,*/ father_index,
                         hostGridPart, subGrid, entities_sharing_same_node, layer, enriched);

              layer += 1;
            }
          }
        }
      }
    }
  } // enrichment

public:
  SubGridList(MacroMicroGridSpecifierType& specifier, bool silent = true)
    : hostSpace_( specifier.fineSpace() )
      , specifier_(specifier)
      , silent_(silent) {
    DSC::Profiler::ScopedTiming st("msfem.subgrid_list");
    DSC_PROFILER.startTiming("msfem.subgrid_list.identify");
    DSC_LOG_INFO << "Starting creation of subgrids." << std::endl << std::endl;


    const HostDiscreteFunctionSpaceType& coarseSpace = specifier_.coarseSpace();

    // the fine grid part
    const HostGridPartType& hostGridPart = hostSpace_.gridPart();

    // the fine grid
    HostGridType& hostGrid = hostSpace_.gridPart().grid();

    const int number_of_nodes = hostGrid.size(2 /*codim*/);    
    
    // -------- identify the entities that share a certain node -------

    std::vector< std::vector< HostEntityPointerType > > entities_sharing_same_node(number_of_nodes);

    // determine the entities that share a common global node with a given index
    for (HostGridEntityIteratorType it = hostSpace_.begin(); it != hostSpace_.end(); ++it)
    {
      int number_of_nodes_in_entity = (*it).template count< 2 >();
      for (int i = 0; i < number_of_nodes_in_entity; i += 1)
      {
        const HostNodePointer node = (*it).template subEntity< 2 >(i);
        int global_index_node = hostGridPart.indexSet().index(*node);

        entities_sharing_same_node[global_index_node].push_back( HostEntityPointerType(*it) );
      }
    }

    // determine the maximum number of oversampling layers
    int max_num_layers = 0;
    for (int i = 0; i < specifier_.getNumOfCoarseEntities(); i += 1)
    {
      if (specifier_.getLayer(i) > max_num_layers)
      { max_num_layers = specifier_.getLayer(i); }
    }

    // the difference in levels between coarse and fine grid
    const int level_difference = specifier_.getLevelDifference();

    // the number of coarse grid entities (of codim 0).
    const int number_of_coarse_grid_entities = specifier_.getNumOfCoarseEntities();

    DSC_LOG_INFO << "number_of_coarse_grid_entities = " << number_of_coarse_grid_entities << std::endl;

    if ( (DSC_CONFIG_GET( "msfem.oversampling_strategy", 1 ) == 2) || (DSC_CONFIG_GET( "msfem.oversampling_strategy", 1 ) == 3) )
    {
      for (int i = 0; i < number_of_coarse_grid_entities; ++i )
      {
         CoarseNodeVectorType coarse_node_vector;
         coarse_node_store_.push_back( coarse_node_vector );
      }
    }      


    //! ----------- create subgrids --------------------

    const HostGridLeafIndexSet& coarseGridLeafIndexSet = coarseSpace.gridPart().grid().leafIndexSet();
    const int oversampling_strategy = DSC_CONFIG_GET( "msfem.oversampling_strategy", 1 );
    EnrichmentMatrixType enriched(boost::extents[number_of_coarse_grid_entities][hostGrid.size(0)][max_num_layers + 1]);
    std::fill( enriched.data(), enriched.data() + enriched.num_elements(), false );
    // loop to initialize subgrids (and to initialize the coarse node vector):
    // -----------------------------------------------------------
    for (HostGridEntityIteratorType coarse_it = coarseSpace.begin(); coarse_it != coarseSpace.end(); ++coarse_it)
    {
      const int coarse_index = coarseGridLeafIndexSet.index(*coarse_it);
      subGridList_.push_back(make_shared<SubGridType>(hostGrid));
      subGridList_[coarse_index]->createBegin();
      
      if ( (oversampling_strategy == 2) || (oversampling_strategy == 3) )
       {
        for (int c = 0; c < coarse_it->geometry().corners(); ++c )
          coarse_node_store_[coarse_index].push_back( coarse_it->geometry().corner(c) );
       }

    }
    // -----------------------------------------------------------
    
    // initialize and fill a vector 'entities_sharing_same_node' that tells you for
    // a given node 'i' which fine grid entities intersect with 'i' 
    // -----------------------------------------------------------
    for (HostGridEntityIteratorType it = hostSpace_.begin(); it != hostSpace_.end(); ++it)
    {
      const int number_of_nodes_in_entity = (*it).template count< 2 >();
      for (int i = 0; i < number_of_nodes_in_entity; i += 1)
      {
        const HostNodePointer node = (*it).template subEntity< 2 >(i);
        const int global_index_node = hostGridPart.indexSet().index(*node);

        entities_sharing_same_node[global_index_node].push_back( HostEntityPointerType(*it) );
      }
    }
    // -----------------------------------------------------------
    
    
    DSC_PROFILER.stopTiming("msfem.subgrid_list.identify");
    DSC_PROFILER.startTiming("msfem.subgrid_list.create");

    // a fine grid iterator for the codim 0 hostgrid entities:
    const HostGridEntityIteratorType host_endit = hostSpace_.end();
    for (HostGridEntityIteratorType host_it = hostSpace_.begin();
         host_it != host_endit;
         ++host_it)
    {
      const HostEntityType& host_entity = *host_it;
      // Dune::Stuff::Grid::printEntity(host_entity);

      // get the coarse-grid-father of host_entity (which is a maxlevel entity)
      HostEntityPointerType level_father_entity = Stuff::Grid::make_father(coarseGridLeafIndexSet,
                                                                           HostEntityPointerType(host_entity),
                                                                           level_difference);
      const int father_index = coarseGridLeafIndexSet.index(*level_father_entity);

      if ( !( subGridList_[father_index]->template contains< 0 >(host_entity) ) )
      { subGridList_[father_index]->insertPartial(host_entity); }
      
      // check the neighbor entities and look if they belong to the same father
      // if yes, continue
      // if not, enrichement with 'n(T)'-layers
      bool all_neighbors_have_same_father = true;
      const HostIntersectionIterator iend = hostGridPart.iend(host_entity);
      for (HostIntersectionIterator iit = hostGridPart.ibegin(host_entity); iit != iend; ++iit)
      {
        if ( iit->neighbor() ) // if there is a neighbor entity
        {
          // check if the neighbor entity is in the subgrid
          const HostEntityPointerType neighborHostEntityPointer = iit->outside();
          HostEntityPointerType level_father_neighbor_entity = Stuff::Grid::make_father(coarseGridLeafIndexSet,
                                                                                        neighborHostEntityPointer,
                                                                                        level_difference);
          if ( !(level_father_neighbor_entity == level_father_entity) )
          {
            all_neighbors_have_same_father = false;
          }
        }
        else {
          all_neighbors_have_same_father = false;
        }
      }

      if (all_neighbors_have_same_father)
      { continue; }

      int layers = specifier_.getLayer(father_index);

      if (layers > 0)
      {
        DSC::Profiler::ScopedTiming st("msfem.subgrid_list.enrichment");
        HostEntityPointerType hep(*host_it);
        enrichment(hep, level_father_entity, father_index, hostGridPart,
                   subGridList_[father_index], entities_sharing_same_node, layers, enriched);
      }
    }


    DSC_PROFILER.startTiming("msfem.subgrid_list.create.finalize");
    for (int i = 0; i < number_of_coarse_grid_entities; ++i)
    {
      assert(int(subGridList_.size()) > i);
      assert(subGridList_[i]);
      subGridList_[i]->createEnd();
      if (!silent_)
      {
        DSC_LOG_INFO << "Subgrid " << i << ":" << std::endl;
        subGridList_[i]->report();
      }

      if (subGridList_[i]->size(2) == 0)
      {
        DSC_LOG_ERROR << "Error." << std::endl
                      << "Error. Created Subgrid with 0 nodes." << std::endl;

        for (HostGridEntityIteratorType coarse_it = specifier_.coarseSpace().begin();
             coarse_it != specifier_.coarseSpace().end(); ++coarse_it)
        {
          const HostEntityType& coarse_entity = *coarse_it;
          const int index = coarseGridLeafIndexSet.index(coarse_entity);
          if (i == index)
          {
            DSC_LOG_ERROR << "We have a problem with the following coarse-grid element:" << std::endl
                          << "coarse element corner(0) = " << coarse_it->geometry().corner(0) << std::endl
                          << "coarse element corner(1) = " << coarse_it->geometry().corner(1) << std::endl
                          << "coarse element corner(2) = " << coarse_it->geometry().corner(2) << std::endl << std::endl;
          }
        }
        DUNE_THROW(Dune::InvalidStateException, "Created Subgrid with 0 nodes");
      }
    }
    DSC_PROFILER.stopTiming("msfem.subgrid_list.create.finalize");
    DSC_PROFILER.stopTiming("msfem.subgrid_list.create");
    //! ----------- end create subgrids --------------------
  }

  SubGridType& getSubGrid(int i) {
    const int size = specifier_.getNumOfCoarseEntities();

    if (i >= size)
    {
      DUNE_THROW(Dune::RangeError, "Error. Subgrid-Index too large.");
    }
    return *(subGridList_[i]);
  } // getSubGrid
  const SubGridType& getSubGrid(int i) const {
    const int size = specifier_.getNumOfCoarseEntities();

    if (i >= size)
    {
      DUNE_THROW(Dune::RangeError, "Error. Subgrid-Index too large.");
    }
    return *(subGridList_[i]);
  } // getSubGrid

  // only required for oversampling strategies with constraints (e.g strategy 2 or 3):
  const CoarseNodeVectorType& getCoarseNodeVector(int i) const {
    const int size = specifier_.getNumOfCoarseEntities();
    if ( DSC_CONFIG_GET( "msfem.oversampling_strategy", 1 ) == 1 )
      DUNE_THROW(Dune::InvalidStateException, "Method 'getCoarseNodeVector' of class 'SubGridList' should not be used in combination with oversampling strategy 1. Check your implementation!");
    if (i >= size)
    {
      DUNE_THROW(Dune::RangeError, "Error. Subgrid-Index too large.");
    }
    return coarse_node_store_[i];
  } // getSubGrid

private:
  const HostDiscreteFunctionSpaceType& hostSpace_;
  MacroMicroGridSpecifierType& specifier_;

  bool silent_;

  typedef std::vector< shared_ptr<SubGridType> > SubGridStorageType;
  SubGridStorageType subGridList_;

  CoarseGridNodeStorageType coarse_node_store_;  

public:
  ~SubGridList()
  {}
};

} //namespace Dune

#endif // #ifndef SUBGRIDLIST_HH
