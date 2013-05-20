#ifdef HAVE_CMAKE_CONFIG
  #include "cmake_config.h"
#elif defined (HAVE_CONFIG_H)
  #include <config.h>
#endif // ifdef HAVE_CMAKE_CONFIG

#include "subgrid-list.hh"

#include <dune/stuff/common/profiler.hh>
#include <dune/stuff/aliases.hh>
#include <dune/stuff/common/logging.hh>
#include <dune/stuff/common/parameter/configcontainer.hh>

#include <dune/multiscale/tools/misc.hh>
#include <dune/multiscale/tools/subgrid_io.hh>

namespace Dune {
namespace Multiscale {
namespace MsFEM {
bool SubGridList::entityPatchInSubgrid(const HostEntityPointerType& hit,
                                          const HostGridPartType& hostGridPart,
                                          shared_ptr< const SubGridType > subGrid,
                                          const EntityPointerCollectionType& entities_sharing_same_node) const
{
  bool patch_in_subgrid = true;

  // loop over the nodes of the enity
  for (int i = 0; i < (*hit).count< 2 >(); ++i) {
    const HostNodePointer node = (*hit).subEntity< 2 >(i);

    const int global_index_node = hostGridPart.indexSet().index(*node);

    for (int j = 0; j < entities_sharing_same_node[global_index_node].size(); ++j) {
      if (!(subGrid->contains< 0 >(*entities_sharing_same_node[global_index_node][j]))) {
        patch_in_subgrid = false;
      }
    }
  }
  return patch_in_subgrid;
} // entityPatchInSubgrid

void SubGridList::enrichment(const HostEntityPointerType& hit,
//                             const HostEntityPointerType& level_father_it,
                             const int& father_index, // father_index = index/number of current subgrid
                             shared_ptr< SubGridType > subGrid,
                             int& layer)
{
  // difference in levels between coarse and fine grid
  const int                      level_difference = specifier_.getLevelDifference();
  HostDiscreteFunctionSpaceType& coarseSpace      = specifier_.coarseSpace();

  const HostGridLeafIndexSet& coarseGridLeafIndexSet = coarseSpace.gridPart().grid().leafIndexSet();

  const HostGridLeafIndexSet& hostGridLeafIndexSet = hostSpace_.gridPart().grid().leafIndexSet();

  for (int l = 0; l <= layer; ++l) {
    enriched_[father_index][hostGridLeafIndexSet.index(*hit)][l] = true;
  }

  //! decrease the number of layers (needed for recursion)
  --layer;

  // loop over the nodes of the fine grid entity
  for (int i = 0; i < (*hit).count< 2 >(); ++i) {
    const HostNodePointer node              = (*hit).subEntity< 2 >(i);
    int                   global_index_node = hostGridPart_.indexSet().index(*node);

    for (size_t j = 0; j < entities_sharing_same_node_[global_index_node].size(); ++j) {
      if (!(subGrid->contains< 0 >(*entities_sharing_same_node_[global_index_node][j]))) {
        subGrid->insertPartial(*entities_sharing_same_node_[global_index_node][j]);
        fine_id_to_subgrid_ids_[hostGridLeafIndexSet.index( *entities_sharing_same_node_[global_index_node][j] )].push_back(father_index);
      
        // get the corners of the father of the fine grid entity 'entities_sharing_same_node_[global_index_node][j]'
        // and add these corners to the vector 'coarse_node_store_[father_index]' (if they are not yet contained)
        if (specifier_.getOversamplingStrategy() == 3) {
          HostEntityPointerType coarse_father = Stuff::Grid::make_father(coarseGridLeafIndexSet,
                                                                         entities_sharing_same_node_[global_index_node][j],
                                                                         level_difference);
          for (int c = 0; c < coarse_father->geometry().corners(); ++c) {
            // ! not an effective search algorithm (should be improved eventually):
            bool node_contained = false;
            for (size_t cn = 0; cn < coarse_node_store_[father_index].size(); ++cn) {
              // hard coding - 2d case:
              if ((coarse_node_store_[father_index][cn][0] == coarse_father->geometry().corner(c)[0])
                  && (coarse_node_store_[father_index][cn][1] == coarse_father->geometry().corner(c)[1])) { 
                node_contained = true; 
              }
            }
            if (!node_contained) { 
              coarse_node_store_[father_index].emplace_back(coarse_father->geometry().corner(c)); 
            }
          }
        }
      }

      if (layer > 0) {
        const int otherEnclosingCoarseCellIndex
                = getEnclosingMacroCellIndex(entities_sharing_same_node_[global_index_node][j]);
        if (father_index!=otherEnclosingCoarseCellIndex) {
          const auto& tmp_entity_ptr = entities_sharing_same_node_[global_index_node][j];
          if (!enriched_[father_index][hostGridLeafIndexSet.index(*tmp_entity_ptr)][layer]) {
            enrichment(tmp_entity_ptr, father_index, subGrid, layer);
            ++layer;
          }
        }
      }
    }
  }
} // enrichment

SubGridList::SubGridList(MacroMicroGridSpecifierType& specifier, bool silent /*= true*/)
  : hostSpace_(specifier.fineSpace()),
    coarseSpace_(specifier.coarseSpace()),
    specifier_(specifier),
    silent_(silent),
    coarseGridLeafIndexSet_(coarseSpace_.gridPart().grid().leafIndexSet()),
    hostGridLeafIndexSet_(hostSpace_.gridPart().grid().leafIndexSet()),
    hostGridPart_(hostSpace_.gridPart()),
    entities_sharing_same_node_(hostGridPart_.grid().size(HostGridPartType::dimension)),
    enriched_(boost::extents[specifier.getNumOfCoarseEntities()][hostGridPart_.grid().size(0)][specifier.maxNumberOverlayLayers() + 1]),
    fineToCoarseMap_(MPIManager::size())

{
  DSC::Profiler::ScopedTiming st("msfem.subgrid_list");

  fine_id_to_subgrid_ids_.resize( hostGridPart_.grid().size(0) );
  
  // initialize the subgrids (no elements are added)
  identifySubGrids();
  
  // add fine grid elements to the subgrids
  createSubGrids();
  
  // finalize the subgrids
  finalizeSubGrids();

}

SubGridList::~SubGridList(){}



/** Get the subgrid belonging to a given coarse cell index.
*
* @param[in] coarseCellIndex The index of a coarse cell.
* @return Returns the subgrid belonging to the coarse cell with the given index.
*/
SubGridList::SubGridType& SubGridList::getSubGrid(int coarseCellIndex)
{
  auto found = subGridList_.find(coarseCellIndex);
  assert(found!=subGridList_.end() && "There is no subgrid for the index you provided!");

  return *(found->second);
} // getSubGrid

/** Get the subgrid belonging to a given coarse cell index.
*
* @param[in] coarseCellIndex The index of a coarse cell.
* @return Returns the subgrid belonging to the coarse cell with the given index.
*/
const SubGridList::SubGridType& SubGridList::getSubGrid(int coarseCellIndex) const
{
  auto found = subGridList_.find(coarseCellIndex);
  assert(found!=subGridList_.end() && "There is no subgrid for the index you provided!");

  return *(found->second);
} // getSubGrid


// given the index of a (codim 0) host grid entity, return the indices of the subgrids that contain the entity
const std::vector< int >& SubGridList::getSubgridIDs_that_contain_entity (int host_enitity_index) const
{
  return fine_id_to_subgrid_ids_[ host_enitity_index ];
}


// only required for oversampling strategies with constraints (e.g strategy 2 or 3):
const SubGridList::CoarseNodeVectorType& SubGridList::getCoarseNodeVector(int i) const
{
  if (specifier_.getOversamplingStrategy() == 1)
    DUNE_THROW(Dune::InvalidStateException,
               "Method 'getCoarseNodeVector' of class 'SubGridList' should not be used in\
                combination with oversampling strategy 1. Check your implementation!");

  if (i >= specifier_.getNumOfCoarseEntities()) {
    DUNE_THROW(Dune::RangeError, "Error. Subgrid-Index too large.");
  }
  return coarse_node_store_[i];
} // getSubGrid

/** Get the index of the coarse cell enclosing the barycentre of a given fine cell.
*
* Given a fine cell, this method computes its barycentre. Using a grid run on the coarse
* grid, it checks which (if any) coarse cell contains the barycentre.
*
* @tparam IteratorType The type of the grid iterator on the coarse grid.
* @param[in] hostEntity The host entity.
* @param[in,out] lastIterator The macro cell that was found in the last run. This should be set to
*                             coarseGrid.begin<0>() in the first run. This iterator will then be
*                             updated and set to the macro element used in this run.
* @param[in] coarseGridLeafIndexSet_ The index set of the coarse grid.
*
*/
int SubGridList::getEnclosingMacroCellIndex(const HostEntityPointerType& hostEntityPointer) {
  // first check, whether we looked for this host entity already
  int myRank = MPIManager::rank();
  int hostEntityIndex = hostGridLeafIndexSet_.index(*hostEntityPointer);
  auto itFound = fineToCoarseMap_[myRank].find(hostEntityIndex);
  if (itFound!=fineToCoarseMap_[myRank].end()) {
    // if so, return the index that was found last time
    return itFound->second;
  }
  static auto lastIterator = coarseSpace_.begin();
  const auto  baryCenter = hostEntityPointer->geometry().center();
  auto macroCellIterator = lastIterator;
  for (; macroCellIterator != coarseSpace_.end(); ++macroCellIterator) {
    if (macroCellIterator->partitionType()==Dune::InteriorEntity) {
      const auto& macroGeo   = macroCellIterator->geometry();
      const auto& refElement = CoarseRefElementType::general(macroGeo.type());

      bool hostEnIsInMacroCell = refElement.checkInside(macroGeo.local(baryCenter));
      if (hostEnIsInMacroCell) {
        lastIterator  = macroCellIterator;
        int macroIndex = coarseGridLeafIndexSet_.index(*macroCellIterator);
        fineToCoarseMap_[myRank][hostEntityIndex] = macroIndex;
        return macroIndex;
      }
    }
  }
  // if we came this far, we did not find the matching enclosing coarse cell for the given
  // fine cell in [lastIterator, coarse grid end]. Start search from beginning
  for (macroCellIterator = coarseSpace_.begin(); macroCellIterator != lastIterator; ++macroCellIterator) {
    if (macroCellIterator->partitionType()==Dune::InteriorEntity) {
      const auto& macroGeo   = macroCellIterator->geometry();
      const auto& refElement = CoarseRefElementType::general(macroGeo.type());
      bool hostEnIsInMacroCell = refElement.checkInside(macroGeo.local(baryCenter));
      if (hostEnIsInMacroCell) {
        lastIterator = macroCellIterator;
        int macroIndex = coarseGridLeafIndexSet_.index(*macroCellIterator);
        fineToCoarseMap_[myRank][hostEntityIndex] = macroIndex;
        return macroIndex;
      }
    }
  }
  // if we came this far, we did not find an enclosing coarse cell at all, issue a warning
  // and return with error code
//  DSC_LOG_DEBUG << "Warning: Host grid entity was not in any coarse grid cell!\n";
  return -1;
}

void SubGridList::identifySubGrids() {
  DSC_PROFILER.startTiming("msfem.subgrid_list.identify");
  DSC_LOG_INFO << "Starting creation of subgrids." << std::endl << std::endl;


  // the fine grid part
  const HostGridPartType& hostGridPart = hostSpace_.gridPart();

  // the fine grid (subgrid needs non-const ref)
  HostGridType& hostGrid = hostSpace_.gridPart().grid();

  // -------- identify the entities that share a certain node -------

  // determine the entities that share a common global node with a given index
  for (const auto& host_entity : hostSpace_) {
    int number_of_nodes_in_entity = host_entity.count< 2 >();
    for (int i = 0; i < number_of_nodes_in_entity; ++i) {
      const HostNodePointer node              = host_entity.subEntity< 2 >(i);
      const int             global_index_node = hostGridPart.indexSet().index(*node);

      entities_sharing_same_node_[global_index_node].emplace_back(host_entity);
    }
  }

  // the number of coarse grid entities (of codim 0).
  const int number_of_coarse_grid_entities = specifier_.getNumOfCoarseEntities();
  const int oversampling_strategy          = specifier_.getOversamplingStrategy();
  DSC_LOG_INFO << "number_of_coarse_grid_entities = " << number_of_coarse_grid_entities << std::endl;

  // determine the maximum number of oversampling layers
  int max_num_layers = specifier_.maxNumberOverlayLayers();

  // the difference in levels between coarse and fine grid
  const int level_difference = specifier_.getLevelDifference();

  if ((oversampling_strategy == 2) || (oversampling_strategy == 3)) {
    coarse_node_store_ = CoarseGridNodeStorageType(number_of_coarse_grid_entities, CoarseNodeVectorType());
  }


  // ! ----------- create subgrids --------------------

  std::fill(enriched_.data(), enriched_.data() + enriched_.num_elements(), false);
  // loop to initialize subgrids (and to initialize the coarse node vector):
  // -----------------------------------------------------------
  for (const auto& coarse_entity : coarseSpace_) {
    // make sure we only create subgrids for interior coarse elements, not
    // for overlap or ghost elements
    assert(coarse_entity.partitionType()==Dune::InteriorEntity);
    const int coarse_index = coarseGridLeafIndexSet_.index(coarse_entity);
    // make sure we did not create a subgrid for the current coarse entity so far
    assert(subGridList_.find(coarse_index)==subGridList_.end());
    subGridList_[coarse_index] = make_shared<SubGridType>(hostGrid);
    subGridList_[coarse_index]->createBegin();

    if ((oversampling_strategy == 2) || (oversampling_strategy == 3)) {
      assert(coarse_index >= 0 && coarse_index < coarse_node_store_.size()
              && "Index set is not suitable for the current implementation!");
      for (int c = 0; c < coarse_entity.geometry().corners(); ++c)
        coarse_node_store_[coarse_index].emplace_back(coarse_entity.geometry().corner(c));
    }
  }

  // -----------------------------------------------------------
  // initialize and fill a vector 'entities_sharing_same_node_' that tells you for
  // a given node 'i' which fine grid entities intersect with 'i'
  // -----------------------------------------------------------
  //! \todo: isn't this exactly the same as in lines 140--149???
  for (HostGridEntityIteratorType it = hostSpace_.begin(); it != hostSpace_.end(); ++it) {
    const auto& localEntity               = *it;
    const int   number_of_nodes_in_entity = localEntity.count< 2 >();
    for (int i = 0; i < number_of_nodes_in_entity; i += 1) {
      const HostNodePointer node              = localEntity.subEntity< 2 >(i);
      const int             global_index_node = hostGridPart.indexSet().index(*node);

      entities_sharing_same_node_[global_index_node].emplace_back(localEntity);
    }
  }
  // -----------------------------------------------------------

  DSC_PROFILER.stopTiming("msfem.subgrid_list.identify");

  return;
}

void SubGridList::createSubGrids() {
  DSC_PROFILER.startTiming("msfem.subgrid_list.create");

  const HostGridLeafIndexSet& hostGridLeafIndexSet = hostSpace_.gridPart().grid().leafIndexSet();
  
  // loop over all host entities and assign them to a macro cell
  auto lastIt = coarseSpace_.begin();
  for (const auto& host_entity : hostSpace_) {
    // get the coarse-grid-father of host_entity (which is a maxlevel entity)...
//    const HostEntityPointerType level_father_entity = Stuff::Grid::make_father(coarseGridLeafIndexSet_,
//            HostEntityPointerType(host_entity),
//            specifier_.getLevelDifference());
    //// ... and its index

    int  macroCellIndex = getEnclosingMacroCellIndex(host_entity);
    // if macroCellIndex is smaller than zero, the enclosing coarse cell was
    // not found. This may be the case if the host cell does not belong to the
    // grid part for the current process.
    if (macroCellIndex>=0) {
      // add host_entity to the subgrid with the index 'macroCellIndex'
      subGridList_[macroCellIndex]->insertPartial(host_entity);
      // add the id of the subgrid to the vecor at position 'index of host grid element'
      fine_id_to_subgrid_ids_[hostGridLeafIndexSet.index( host_entity )].push_back(macroCellIndex);
      
      // check the neighbor entities and look if they belong to the same father
      // if yes, continue
      // if not, enrichment with 'n(T)' layers
      bool                           all_neighbors_have_same_father = true;
      const HostIntersectionIterator iend                           = hostGridPart_.iend(host_entity);
      for (HostIntersectionIterator iit = hostGridPart_.ibegin(host_entity); (iit != iend) && all_neighbors_have_same_father; ++iit) {
        if (iit->neighbor()) {
          // if there is a neighbor entity
          // check if the neighbor entity is in the subgrid
          if (getEnclosingMacroCellIndex(iit->outside()) != macroCellIndex)
            all_neighbors_have_same_father = false;
        } else {
          all_neighbors_have_same_father = false;
        }
      }

      if (!all_neighbors_have_same_father) {
        int layers = specifier_.getNoOfLayers(macroCellIndex);
        if (layers > 0) {
          DSC::Profiler::ScopedTiming enrichment_st("msfem.subgrid_list.enrichment");
          const HostEntityPointerType hep(host_entity);
          enrichment(hep, macroCellIndex, subGridList_[macroCellIndex], layers);
        }
      }
    }
  }
  DSC_PROFILER.stopTiming("msfem.subgrid_list.create");

  return;
}

void SubGridList::finalizeSubGrids() {
  DSC_PROFILER.startTiming("msfem.subgrid_list.create.finalize");
  int i=0;
  for (auto subGridIt : subGridList_) {
    // finish the creation of each subgrid
    subGridIt.second->createEnd();

    // report some infos about the subgrids if desired
    if (!silent_) {
      DSC_LOG_INFO << "Subgrid " << i << ":" << std::endl;
      subGridIt.second->report();
    }

    // error handling
    if (subGridIt.second->size(2) == 0) {
      DSC_LOG_ERROR << "Error." << std::endl
              << "Error. Created Subgrid with 0 nodes." << std::endl;

      for (const auto& coarse_it : coarseSpace_) {
        const int             index         = coarseGridLeafIndexSet_.index(coarse_it);
        if (subGridIt.first == index) {
          DSC_LOG_ERROR << "We have a problem with the following coarse-grid element:" << std::endl
                  << "coarse element corner(0) = " << coarse_it.geometry().corner(0) << std::endl
                  << "coarse element corner(1) = " << coarse_it.geometry().corner(1) << std::endl
                  << "coarse element corner(2) = " << coarse_it.geometry().corner(2) << std::endl << std::endl;
        }
      }
      DUNE_THROW(Dune::InvalidStateException, "Created Subgrid with 0 nodes");
    }
  }
  DSC_PROFILER.stopTiming("msfem.subgrid_list.create.finalize");

  return;
}

} // namespace MsFEM {
} // namespace Multiscale {
} // namespace Dune {