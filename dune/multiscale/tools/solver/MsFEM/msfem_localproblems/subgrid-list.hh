#ifndef SUBGRIDLIST_HH
#define SUBGRIDLIST_HH

#include <boost/noncopyable.hpp>

#include <dune/common/fmatrix.hh>
#include <dune/common/exceptions.hh>
#include <dune/subgrid/subgrid.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/operator/2order/lagrangematrixsetup.hh>
#include <dune/multiscale/tools/misc.hh>

// / done

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

  //! ---------------- typedefs for the SubgridDiscreteFunctionSpace -----------------------
  // ( typedefs for the local grid and the corresponding local ('sub') )discrete space )

  //! type of grid
  typedef SubGridImp SubGridType;

  //! type of grid part
  typedef LeafGridPart< SubGridType > SubGridPartType;

private:
  template< typename EntityPointerCollectionType >
  bool entity_patch_in_subgrid(const HostEntityPointerType& hit,
                               const HostGridPartType& hostGridPart,
                               const SubGridType& subGrid,
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
                  const int& father_index,
                  const HostGridPartType& hostGridPart,
                  SubGridType& subGrid,
                  EntityPointerCollectionType& entities_sharing_same_node,
                  int& layer,
                  bool***& enriched) {
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

    // loop over the nodes of the enity
    for (int i = 0; i < (*hit).template count< 2 >(); i += 1)
    {
      const HostNodePointer node = (*hit).template subEntity< 2 >(i);
      int global_index_node = hostGridPart.indexSet().index(*node);

      for (size_t j = 0; j < entities_sharing_same_node[global_index_node].size(); j += 1)
      {
        if ( !( subGrid.template contains< 0 >(*entities_sharing_same_node[global_index_node][j]) ) )
        {
          subGrid.insertPartial(*entities_sharing_same_node[global_index_node][j]);
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
    DSC_LOG_INFO << "Starting creation of subgrids." << std::endl << std::endl;

    const HostDiscreteFunctionSpaceType& coarseSpace = specifier_.coarseSpace();

    const HostGridPartType& hostGridPart = hostSpace_.gridPart();

    HostGridType& hostGrid = hostSpace_.gridPart().grid();

    const int number_of_nodes = hostGrid.size(2 /*codim*/);

    // -------- identify the entities that share a certain node -------

    std::vector< std::vector< HostEntityPointerType > > entities_sharing_same_node(number_of_nodes);

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

    int max_num_layers = 0;
    for (int i = 0; i < specifier_.getNumOfCoarseEntities(); i += 1)
    {
      if (specifier_.getLayer(i) > max_num_layers)
      { max_num_layers = specifier_.getLayer(i); }
    }

    // difference in levels between coarse and fine grid
    const int level_difference = specifier_.getLevelDifference();

    // number of coarse grid entities (of codim 0).
    const int number_of_coarse_grid_entities = specifier_.getNumOfCoarseEntities();

    DSC_LOG_INFO << "number_of_coarse_grid_entities = " << number_of_coarse_grid_entities << std::endl;

    subGridList_ = new SubGridType *[number_of_coarse_grid_entities];

    //! ----------- create subgrids --------------------

    const HostGridLeafIndexSet& coarseGridLeafIndexSet = coarseSpace.gridPart().grid().leafIndexSet();

    for (HostGridEntityIteratorType coarse_it = coarseSpace.begin(); coarse_it != coarseSpace.end(); ++coarse_it)
    {
      const int coarse_index = coarseGridLeafIndexSet.index(*coarse_it);

      subGridList_[coarse_index] = new SubGridType(hostGrid);
      subGridList_[coarse_index]->createBegin();
    }

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

    bool*** enriched = new bool**[number_of_coarse_grid_entities];
    for (int k = 0; k < number_of_coarse_grid_entities; k += 1)
    {
      enriched[k] = new bool*[hostGrid.size(0 /*codim*/)];
      for (int m = 0; m < hostGrid.size(0 /*codim*/); m += 1)
      {
        enriched[k][m] = new bool[max_num_layers + 1];
      }
    }

    for (int k = 0; k < number_of_coarse_grid_entities; k += 1)
    {
      for (int m = 0; m < hostGrid.size(0 /*codim*/); m += 1)
      {
        for (int l = 0; l < max_num_layers + 1; l += 1)
          enriched[k][m][l] = false;
      }
    }

    // a fine grid iterator for the codim 0 hostgrid entities:
    const HostGridEntityIteratorType host_endit = hostSpace_.end();
    for (HostGridEntityIteratorType host_it = hostSpace_.begin();
         host_it != host_endit;
         ++host_it)
    {
      const HostEntityType& host_entity = *host_it;

      const int DUNE_UNUSED(number_of_nodes_in_entity) = (*host_it).template count< 2 >();

      // get the coarse-grid-father of host_entity (which is a maxlevel entity)
      HostEntityPointerType level_father_entity = Stuff::Grid::make_father(coarseGridLeafIndexSet,
                                                                           HostEntityPointerType(*host_it),
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
          const HostEntityType& DUNE_UNUSED(neighborHostEntity) = *neighborHostEntityPointer;

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
        HostEntityPointerType hep(*host_it);
        enrichment(hep, level_father_entity, father_index, hostGridPart,
                   *subGridList_[father_index], entities_sharing_same_node, layers, enriched);
      }
    }

    for (int i = 0; i < number_of_coarse_grid_entities; ++i)
    {
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

private:
  const HostDiscreteFunctionSpaceType& hostSpace_;
  MacroMicroGridSpecifierType& specifier_;

  bool silent_;

  SubGridType** subGridList_;

public:
  // Destruktor (Dynamisch angeforderter Speicher wird wieder freigegeben)
  ~SubGridList() {
    const int size = specifier_.getNumOfCoarseEntities();

    for (int i = 0; i < size; ++i)
      delete subGridList_[i];
    delete[] subGridList_;
  }
};
}

#endif // #ifndef SUBGRIDLIST_HH
