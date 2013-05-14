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
                             const HostEntityPointerType& level_father_it,
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
        const HostEntityPointerType father = Stuff::Grid::make_father(coarseGridLeafIndexSet,
                                                                      entities_sharing_same_node_[global_index_node][j],
                                                                      level_difference);
        if (father != level_father_it) {
          const auto& tmp_entity_ptr = entities_sharing_same_node_[global_index_node][j];
          if (!enriched_[father_index][hostGridLeafIndexSet.index(*tmp_entity_ptr)][layer]) {
            enrichment(tmp_entity_ptr, level_father_it, father_index, subGrid, layer);

            layer += 1;
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
    hostGridPart_(hostSpace_.gridPart()),
    entities_sharing_same_node_(hostGridPart_.grid().size(HostGridPartType::dimension)),
    enriched_(boost::extents[specifier.getNumOfCoarseEntities()][hostGridPart_.grid().size(0)][specifier.maxNumberOverlayLayers() + 1])

{
  DSC::Profiler::ScopedTiming st("msfem.subgrid_list");

  identifySubGrids();
  createSubGrids();
  finalizeSubGrids();
}

SubGridList::~SubGridList(){}



SubGridList::SubGridType& SubGridList::getSubGrid(int i)
{
  const int size = specifier_.getNumOfCoarseEntities();

  if (i >= size) {
    DUNE_THROW(Dune::RangeError, "Error. Subgrid-Index too large.");
  }
  return *(subGridList_[i]);
} // getSubGrid

const SubGridList::SubGridType& SubGridList::getSubGrid(int i) const
{
  const int size = specifier_.getNumOfCoarseEntities();

  if (i >= size) {
    DUNE_THROW(Dune::RangeError, "Error. Subgrid-Index too large.");
  }
  return *(subGridList_[i]);
} // getSubGrid

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
template<class IteratorType>
int SubGridList::getEnclosingMacroCellIndex(const HostEntityPointerType& hostEntityPointer,
        IteratorType& lastIterator) {
  const auto  baryCenter = hostEntityPointer->geometry().center();
  IteratorType macroCellIterator = lastIterator;
  for (; macroCellIterator != coarseSpace_.end(); ++macroCellIterator) {
    const auto& macroGeo   = macroCellIterator->geometry();
    const auto& refElement = CoarseRefElementType::general(macroGeo.type());

    bool hostEnIsInMacroCell = refElement.checkInside(macroGeo.local(baryCenter));
    if (hostEnIsInMacroCell) {
      lastIterator  = macroCellIterator;
      return coarseGridLeafIndexSet_.index(*macroCellIterator);
    }
  }
  // if we came this far, we did not find the matching enclosing coarse cell for the given
  // fine cell in [lastIterator, coarse grid end]. Start search from beginning
  for (macroCellIterator = coarseSpace_.begin(); macroCellIterator != lastIterator; ++macroCellIterator) {
    const auto& macroGeo   = macroCellIterator->geometry();
    const auto& refElement = CoarseRefElementType::general(macroGeo.type());
    bool hostEnIsInMacroCell = refElement.checkInside(macroGeo.local(baryCenter));
    if (hostEnIsInMacroCell) {
      lastIterator = macroCellIterator;
      return coarseGridLeafIndexSet_.index(*macroCellIterator);
    }
  }
  // if we came this far, we did not find an enclosing coarse cell at all, issue a warning
  // and return with error code
  DSC_LOG_INFO << "Warning: Host grid entity was not in any coarse grid cell!\n";
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
    const int coarse_index = coarseGridLeafIndexSet_.index(coarse_entity);
    subGridList_.emplace_back(new SubGridType(hostGrid));
    subGridList_.back()->createBegin();

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

  // loop over all host entities and assign them to a macro cell
  auto lastIt = coarseSpace_.begin();
  for (const auto& host_entity : hostSpace_) {
    // get the coarse-grid-father of host_entity (which is a maxlevel entity)...
    const HostEntityPointerType level_father_entity = Stuff::Grid::make_father(coarseGridLeafIndexSet_,
            HostEntityPointerType(host_entity),
            specifier_.getLevelDifference());
    //// ... and its index
    int  macroCellIndex = getEnclosingMacroCellIndex(host_entity, lastIt);
    subGridList_[macroCellIndex]->insertPartial(host_entity);

    // check the neighbor entities and look if they belong to the same father
    // if yes, continue
    // if not, enrichment with 'n(T)' layers
    bool                           all_neighbors_have_same_father = true;
    const HostIntersectionIterator iend                           = hostGridPart_.iend(host_entity);
    for (HostIntersectionIterator iit = hostGridPart_.ibegin(host_entity); (iit != iend) && all_neighbors_have_same_father; ++iit) {
      if (iit->neighbor()) {
        // if there is a neighbor entity
        // check if the neighbor entity is in the subgrid
        int neighbourEnlosingMacroCellIndex = getEnclosingMacroCellIndex(iit->outside(), lastIt);
        if (neighbourEnlosingMacroCellIndex!=macroCellIndex)
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
        enrichment(hep, level_father_entity, macroCellIndex, subGridList_[macroCellIndex], layers);
      }
    }
  }
  DSC_PROFILER.stopTiming("msfem.subgrid_list.create");

  return;
}

void SubGridList::finalizeSubGrids() {
  DSC_PROFILER.startTiming("msfem.subgrid_list.create.finalize");
  const int numCoarseEntities = specifier_.getNumOfCoarseEntities();
  for (int i = 0; i < numCoarseEntities; ++i) {
    assert(int(subGridList_.size()) > i);
    assert(subGridList_[i]);

    // finish the creation of each subgrid
    subGridList_[i]->createEnd();

    // report some infos about the subgrids if desired
    if (!silent_) {
      DSC_LOG_INFO << "Subgrid " << i << ":" << std::endl;
      subGridList_[i]->report();
    }

    // error handling
    if (subGridList_[i]->size(2) == 0) {
      DSC_LOG_ERROR << "Error." << std::endl
              << "Error. Created Subgrid with 0 nodes." << std::endl;

      for (HostGridEntityIteratorType coarse_it = specifier_.coarseSpace().begin();
           coarse_it != specifier_.coarseSpace().end(); ++coarse_it) {
        const HostEntityType& coarse_entity = *coarse_it;
        const int             index         = coarseGridLeafIndexSet_.index(coarse_entity);
        if (i == index) {
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

  return;
}

} // namespace MsFEM {
} // namespace Multiscale {
} // namespace Dune {