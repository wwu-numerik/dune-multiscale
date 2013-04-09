// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef SUBGRIDLIST_HH
#define SUBGRIDLIST_HH

#include <boost/noncopyable.hpp>
#include <boost/multi_array.hpp>

#include <dune/common/fmatrix.hh>
#include <dune/common/shared_ptr.hh>
#include <dune/common/exceptions.hh>
#include <dune/stuff/grid/entity.hh>

#include <dune/multiscale/tools/subgrid_io.hh>
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
namespace Multiscale {
namespace MsFEM {

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
  //! type of grid partition
  typedef typename HostDiscreteFunctionSpaceType::GridPartType HostGridPartType;

  //! type of grid
private:
  typedef typename HostDiscreteFunctionSpaceType::GridType HostGridType;
  typedef typename HostGridType::Traits::LeafIndexSet HostGridLeafIndexSet;
  typedef typename HostDiscreteFunctionSpaceType::IteratorType HostGridEntityIteratorType;
  typedef typename HostGridEntityIteratorType::Entity HostEntityType;
  typedef typename HostEntityType::EntityPointer HostEntityPointerType;
  typedef typename HostEntityType::template Codim< 2 >::EntityPointer HostNodePointer;
  typedef typename HostGridPartType::IntersectionIteratorType HostIntersectionIterator;

  //! type of (non-discrete )function space
  typedef typename HostDiscreteFunctionSpaceType::FunctionSpaceType FunctionSpaceType;

  //! type of domain
  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef std::vector< DomainType > CoarseNodeVectorType;
  typedef std::vector< CoarseNodeVectorType > CoarseGridNodeStorageType;
  typedef boost::multi_array<bool, 3> EnrichmentMatrixType;

public:
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
          if ( specifier_.getOversamplingStrategy() == 3 )
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
              { coarse_node_store_[father_index].emplace_back( coarse_father->geometry().corner(c) ); }
            }
          }
        }

        if (layer > 0)
        {
          const HostEntityPointerType father = Stuff::Grid::make_father(coarseGridLeafIndexSet,
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

    // the fine grid (subgrid needs non-const ref
    HostGridType& hostGrid = hostSpace_.gridPart().grid();

    const int number_of_nodes = hostGrid.size(2 /*codim*/);    
    
    // -------- identify the entities that share a certain node -------

    std::vector< std::vector< HostEntityPointerType > > entities_sharing_same_node(number_of_nodes);

    // determine the entities that share a common global node with a given index
    for (const auto& host_entity : hostSpace_)
    {
      int number_of_nodes_in_entity = host_entity.template count< 2 >();
      for (int i = 0; i < number_of_nodes_in_entity; i += 1)
      {
        const HostNodePointer node = host_entity.template subEntity< 2 >(i);
        const int global_index_node = hostGridPart.indexSet().index(*node);

        entities_sharing_same_node[global_index_node].emplace_back(host_entity);
      }
    }

    // determine the maximum number of oversampling layers
    int max_num_layers = 0;
    for (int i = 0; i < specifier_.getNumOfCoarseEntities(); i += 1)
    {
      max_num_layers = std::max(max_num_layers, specifier_.getLayer(i));
    }

    // the difference in levels between coarse and fine grid
    const int level_difference = specifier_.getLevelDifference();

    // the number of coarse grid entities (of codim 0).
    const int number_of_coarse_grid_entities = specifier_.getNumOfCoarseEntities();
    const int oversampling_strategy = specifier_.getOversamplingStrategy();
    DSC_LOG_INFO << "number_of_coarse_grid_entities = " << number_of_coarse_grid_entities << std::endl;

    if ( (oversampling_strategy == 2) || (oversampling_strategy == 3) )
    {
      coarse_node_store_ = CoarseGridNodeStorageType(number_of_coarse_grid_entities, CoarseNodeVectorType());
    }      


    //! ----------- create subgrids --------------------

    const HostGridLeafIndexSet& coarseGridLeafIndexSet = coarseSpace.gridPart().grid().leafIndexSet();
    EnrichmentMatrixType enriched(boost::extents[number_of_coarse_grid_entities][hostGrid.size(0)][max_num_layers + 1]);
    std::fill( enriched.data(), enriched.data() + enriched.num_elements(), false );
    // loop to initialize subgrids (and to initialize the coarse node vector):
    // -----------------------------------------------------------
    for (const auto& coarse_entity : coarseSpace)
    {
      const int coarse_index = coarseGridLeafIndexSet.index(coarse_entity);
      subGridList_.emplace_back(new SubGridType(hostGrid));
      subGridList_[coarse_index]->createBegin();
      
      if ( (oversampling_strategy == 2) || (oversampling_strategy == 3) )
       {
        for (int c = 0; c < coarse_entity.geometry().corners(); ++c )
          coarse_node_store_[coarse_index].emplace_back( coarse_entity.geometry().corner(c) );
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

        entities_sharing_same_node[global_index_node].emplace_back(*it);
      }
    }
    // -----------------------------------------------------------
    
    
    DSC_PROFILER.stopTiming("msfem.subgrid_list.identify");
    DSC_PROFILER.startTiming("msfem.subgrid_list.create");

    for (const auto& host_entity : hostSpace_)
    {
      // Dune::Stuff::Grid::printEntity(host_entity);

      // get the coarse-grid-father of host_entity (which is a maxlevel entity)
      const HostEntityPointerType level_father_entity = Stuff::Grid::make_father(coarseGridLeafIndexSet,
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
          const HostEntityPointerType level_father_neighbor_entity = Stuff::Grid::make_father(coarseGridLeafIndexSet,
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
        DSC::Profiler::ScopedTiming enrichment_st("msfem.subgrid_list.enrichment");
        const HostEntityPointerType hep(host_entity);
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
    if ( specifier_.getOversamplingStrategy() == 1 )
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

} //namespace MsFEM {
} //namespace Multiscale {
} //namespace Dune {

#endif // #ifndef SUBGRIDLIST_HH
