#include <config.h>
#include <assert.h>
#include <parallel/algorithm>
#include <boost/assert.hpp>
#include <boost/multi_array/multi_array_ref.hpp>
#include <dune/common/exceptions.hh>
#include <dune/fem/misc/mpimanager.hh>
#include <dune/stuff/common/logging.hh>
#include <dune/stuff/common/profiler.hh>
#include <dune/stuff/common/ranges.hh>
#include <dune/stuff/common/float_cmp.hh>
#include <algorithm>
#include <iterator>
#include <ostream>
#include <utility>
#include <memory>
#include <dune/multiscale/tools/misc.hh>

#include "subgrid-list.hh"

namespace Dune {
namespace Multiscale {
namespace MsFEM {
bool SubGridList::entityPatchInSubgrid(const HostEntityPointerType& hit, const HostGridPartType& hostGridPart,
                                       shared_ptr<const SubGridType> subGrid,
                                       const EntityPointerCollectionType& entities_sharing_same_node) const {
  bool patch_in_subgrid = true;

  // loop over the nodes of the enity
  for (int i = 0; i < (*hit).count<HostGridType::dimension>(); ++i) {
    const HostNodePointer node = (*hit).subEntity<HostGridType::dimension>(i);

    const int global_index_node = hostGridPart.indexSet().index(*node);

    for (std::size_t j = 0; j < entities_sharing_same_node[global_index_node].size(); ++j) {
      if (!(subGrid->contains<0>(*entities_sharing_same_node[global_index_node][j]))) {
        patch_in_subgrid = false;
      }
    }
  }
  return patch_in_subgrid;
} // entityPatchInSubgrid

void SubGridList::enrichment(const HostEntityPointerType& hit,
                             //                             const HostEntityPointerType& level_father_it,
                             const std::size_t& subgrid_index, // subgrid_index = father_index = index/number of current
                                                               // subgrid
                             std::shared_ptr<SubGridType> subGrid, std::size_t& layer) {
  assert(layer >= 1);
  // difference in levels between coarse and fine grid
  const auto level_difference = specifier_.getLevelDifference();
  const auto& coarseSpace = specifier_.coarseSpace();
  const auto& coarseGridLeafIndexSet = coarseSpace.gridPart().grid().leafIndexSet();
  const auto& hostGridLeafIndexSet = hostSpace_.gridPart().grid().leafIndexSet();

  for (auto l : DSC::valueRange(layer)) {
    enriched_[subgrid_index][hostGridLeafIndexSet.index(*hit)][l] = true;
  }

  //! decrease the number of layers (needed for recursion)
  --layer;

  // loop over the nodes of the fine grid entity
  for (int i = 0; i < (*hit).count<HostGridType::dimension>(); ++i) {
    const HostNodePointer node = (*hit).subEntity<HostGridType::dimension>(i);
    int global_index_node = hostGridPart_.indexSet().index(*node);

    // loop over the the fine grid entities that share the node
    for (size_t j = 0; j < entities_sharing_same_node_[global_index_node].size(); ++j) {
      // if the subgrid does not yet contain the fine grid entity ..
      if (!(subGrid->contains<0>(*entities_sharing_same_node_[global_index_node][j]))) {
        // .. add it to the subgrid
        subGrid->insertPartial(*entities_sharing_same_node_[global_index_node][j]);
        // also add the information that the fine grid element is contained in the subgrid
        fine_id_to_subgrid_ids_[hostGridLeafIndexSet.index(*entities_sharing_same_node_[global_index_node][j])]
            .push_back(subgrid_index);

        // get the corners of the father of the fine grid entity 'entities_sharing_same_node_[global_index_node][j]'
        // and add these corners to the vector 'coarse_node_store_[subgrid_index]' (if they are not yet contained)
        if (specifier_.getOversamplingStrategy() == 3) {

          HostEntityPointerType& current_fine_entity = entities_sharing_same_node_[global_index_node][j];
          HostEntityPointerType coarse_father =
              DSG::make_father(coarseGridLeafIndexSet, current_fine_entity, level_difference);
          for (int c = 0; c < coarse_father->geometry().corners(); ++c) {

            // add the corners of the coarse_father to the 'extended_coarse_node_store_' if they are not yet contained
            bool coarse_node_contained = false;
            // check if node is already contained in the vector 'extended_coarse_node_store_[subgrid_index]'
            for (size_t cn = 0; cn < extended_coarse_node_store_[subgrid_index].size(); ++cn) {
              // ! not an effective search algorithm (should be improved eventually):
              for (unsigned int d = 0; d < DomainType::dimension; ++d) {
                if (DSC::FloatCmp::ne(extended_coarse_node_store_[subgrid_index][cn][d],
                                      coarse_father->geometry().corner(c)[d])) {
                  break;
                } else if (d == (DomainType::dimension - 1)) {
                  coarse_node_contained = true;
                }
              }
            }
            if (!coarse_node_contained) {
              extended_coarse_node_store_[subgrid_index].emplace_back(coarse_father->geometry().corner(c));
            }

            for (int c_fine = 0; c_fine < current_fine_entity->geometry().corners(); ++c_fine) {

              // check if one of the corners of the fine entity is identical to one of the corners of the coarse father
              // if we find such a corner, (and if it is not already contained) it must be added to the
              // coarse_node_store_
              bool fine_corner_is_coarse_corner = false;
              for (unsigned int d = 0; d < DomainType::dimension; ++d) {
                if (DSC::FloatCmp::ne(current_fine_entity->geometry().corner(c_fine)[d],
                                      coarse_father->geometry().corner(c)[d])) {
                  break;
                } else if (d == (DomainType::dimension - 1)) {
                  fine_corner_is_coarse_corner = true;
                }
              }

              if (fine_corner_is_coarse_corner) {

                bool node_contained =
                    false; // check if node is already contained in the vector 'coarse_node_store_[subgrid_index]'
                for (size_t cn = 0; cn < coarse_node_store_[subgrid_index].size(); ++cn) {

                  for (unsigned int d = 0; d < DomainType::dimension; ++d) {
                    if (DSC::FloatCmp::ne(coarse_node_store_[subgrid_index][cn][d],
                                          coarse_father->geometry().corner(c)[d])) {
                      break;
                    } else if (d == (DomainType::dimension - 1)) {
                      node_contained = true;
                    }
                  }
                }
                if (!node_contained) {
                  coarse_node_store_[subgrid_index].emplace_back(coarse_father->geometry().corner(c));
                }
              }
            }
          }
        }
      }

      if (layer > 0) {
        const auto otherEnclosingCoarseCellIndex =
            getEnclosingMacroCellIndex(entities_sharing_same_node_[global_index_node][j]);
        if (subgrid_index != otherEnclosingCoarseCellIndex) {
          const auto& tmp_entity_ptr = entities_sharing_same_node_[global_index_node][j];
          if (!enriched_[subgrid_index][hostGridLeafIndexSet.index(*tmp_entity_ptr)][layer]) {
            enrichment(tmp_entity_ptr, subgrid_index, subGrid, layer);
            ++layer;
          }
        }
      }
    }
  }
} // enrichment

SubGridList::SubGridList(MacroMicroGridSpecifierType& specifier, bool silent /*= true*/)
  : hostSpace_(specifier.fineSpace())
  , coarseSpace_(specifier.coarseSpace())
  , specifier_(specifier)
  , silent_(silent)
  , coarseGridLeafIndexSet_(coarseSpace_.gridPart().grid().leafIndexSet())
  , hostGridLeafIndexSet_(hostSpace_.gridPart().grid().leafIndexSet())
  , hostGridPart_(hostSpace_.gridPart())
  , entities_sharing_same_node_(hostGridPart_.grid().size(HostGridPartType::dimension))
  , enriched_(boost::extents[specifier.getNumOfCoarseEntities()][hostGridPart_.grid().size(0)]
                            [specifier.maxNumberOverlayLayers() + 1])
  , fineToCoarseMap_(Fem::MPIManager::size())
  , threadIterators_(coarseSpace_.gridPart())
{
  DSC::Profiler::ScopedTiming st("msfem.subgrid_list");

  threadIterators_.update();
  fine_id_to_subgrid_ids_.resize(hostGridPart_.grid().size(0));

  //! @todo temp!
  for (const auto& hostEntity : DSC::viewRange(hostGridPart_.grid().leafView())) {
    getEnclosingMacroCellId(hostEntity);
  }

  // initialize the subgrids (no elements are added)
  identifySubGrids();

  // add fine grid elements to the subgrids
  createSubGrids();

  // finalize the subgrids
  finalizeSubGrids();
}

SubGridList::~SubGridList() {}

/** Get the subgrid belonging to a given coarse cell index.
*
* @param[in] coarseCellIndex The index of a coarse cell.
* @return Returns the subgrid belonging to the coarse cell with the given index.
*/
SubGridList::SubGridType& SubGridList::getSubGrid(std::size_t coarseCellIndex) {
  auto found = subGridList_.find(coarseCellIndex);
  BOOST_ASSERT_MSG(found != subGridList_.end(), "There is no subgrid for the index you provided!");
  assert(found->second);
  return *(found->second);
} // getSubGrid

/** Get the subgrid belonging to a given coarse cell index.
*
* @param[in] coarseCellIndex The index of a coarse cell.
* @return Returns the subgrid belonging to the coarse cell with the given index.
*/
const SubGridList::SubGridType& SubGridList::getSubGrid(std::size_t coarseCellIndex) const {
  auto found = subGridList_.find(coarseCellIndex);
  BOOST_ASSERT_MSG(found != subGridList_.end(), "There is no subgrid for the index you provided!");
  assert(found->second);
  return *(found->second);
} // getSubGrid

/** Get the subgrid belonging to a given coarse cell.
*
* @param[in] coarseCell The coarse cell.
* @return Returns the subgrid belonging to the given coarse cell.
*/
const SubGridList::SubGridType& SubGridList::getSubGrid(const CoarseEntityType& entity) const {
  const int index = coarseGridLeafIndexSet_.index(entity);
  return getSubGrid(index);
} // getSubGrid

/** Get the subgrid belonging to a given coarse cell.
*
* @param[in] coarseCell The coarse cell.
* @return Returns the subgrid belonging to the given coarse cell.
*/
SubGridList::SubGridType& SubGridList::getSubGrid(const CoarseEntityType& entity) {
  const int index = coarseGridLeafIndexSet_.index(entity);
  return getSubGrid(index);
} // getSubGrid

const SubGridList::EntityPointerCollectionType& SubGridList::getNodeEntityMap() { return entities_sharing_same_node_; }

// given the id of a subgrid, return the entity seed for the 'base coarse entity'
// (i.e. the coarse entity that the subgrid was constructed from by enrichment )
const SubGridList::CoarseGridEntitySeed& SubGridList::get_coarse_entity_seed(std::size_t i) const {
  // the following returns the mapped element for index i if present,
  // if not, an out-of-range exception is thrown
  assert(false);//need to eliminate narrowing conversion
  return subgrid_id_to_base_coarse_entity_.at(i);
}

// given the index of a (codim 0) host grid entity, return the indices of the subgrids that contain the entity
const std::vector<std::size_t>& SubGridList::getSubgridIDs_that_contain_entity(std::size_t host_enitity_index) const {
  return fine_id_to_subgrid_ids_[host_enitity_index];
}

// only required for oversampling strategies with constraints (e.g strategy 2 or 3):
// for each given subgrid index return the vecor of ALL coarse nodes (global coordinates) that are in the subgrid,
// this also includes the coarse nodes on the boundary of U(T), even if this is a global Dirichlet node!
const SubGridList::CoarseNodeVectorType& SubGridList::getCoarseNodeVector(std::size_t i) const {
  if (specifier_.getOversamplingStrategy() == 1)
    DUNE_THROW(Dune::InvalidStateException, "Method 'getCoarseNodeVector' of class 'SubGridList' should not be used in\
                combination with oversampling strategy 1. Check your implementation!");

  if (i >= specifier_.getNumOfCoarseEntities()) {
    DUNE_THROW(Dune::RangeError, "Error. Subgrid-Index too large.");
  }
  return coarse_node_store_[i];
} // getSubGrid

// only required for oversampling strategy 3:
// this method only differs from the method 'getCoarseNodeVector' if the oversampling patch
// cannot be excluively described by a union of coarse grid elements.
// According to the definition of the LOD 'not full coarse layers' require that the averaging
// property of the weighted Clement operator is also applied to those coarse nodes, where
// the corresponding basis function has a nonempty intersection with the patch
const SubGridList::CoarseNodeVectorType& SubGridList::getExtendedCoarseNodeVector(std::size_t i) const {
  if ((specifier_.getOversamplingStrategy() == 1) || (specifier_.getOversamplingStrategy() == 2))
    DUNE_THROW(Dune::InvalidStateException,
               "Method 'getExtendendCoarseNodeVector' of class 'SubGridList' should not be used in\
                combination with oversampling strategy 1 or 2. Check your implementation!");

  if (i >= specifier_.getNumOfCoarseEntities()) {
    DUNE_THROW(Dune::RangeError, "Error. Subgrid-Index too large.");
  }
  return extended_coarse_node_store_[i];
} // getSubGrid

// get number of sub grids
std::size_t SubGridList::getNumberOfSubGrids() const { return specifier_.getNumOfCoarseEntities(); }
std::size_t SubGridList::size() const { return specifier_.getNumOfCoarseEntities(); }

SubGridList::SubGridPartType SubGridList::gridPart(std::size_t i) { return SubGridPartType(getSubGrid(i)); }

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
std::size_t SubGridList::getEnclosingMacroCellIndex(const HostEntityPointerType& hostEntityPointer) {
  // first check, whether we looked for this host entity already
  const auto myRank = Fem::MPIManager::rank();
  const auto hostEntityIndex = hostGridLeafIndexSet_.index(*hostEntityPointer);
  const auto itFound = fineToCoarseMap_[myRank].find(hostEntityIndex);
  if (itFound != fineToCoarseMap_[myRank].end()) {
    // if so, return the index that was found last time
    return itFound->second;
  }
  static auto lastIterator = coarseSpace_.begin();
  const auto baryCenter = hostEntityPointer->geometry().center();
  auto macroCellIterator = lastIterator;
  for (; macroCellIterator != coarseSpace_.end(); ++macroCellIterator) {
    if (macroCellIterator->partitionType() == Dune::InteriorEntity) {
      const auto& macroGeo = macroCellIterator->geometry();
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
  // if we came this far, we did not find the matching enclosing coarse cell for the given
  // fine cell in [lastIterator, coarse grid end]. Start search from beginning
  for (macroCellIterator = coarseSpace_.begin(); macroCellIterator != lastIterator; ++macroCellIterator) {
    if (macroCellIterator->partitionType() == Dune::InteriorEntity) {
      const auto& macroGeo = macroCellIterator->geometry();
      const auto& refElement = CoarseRefElementType::general(macroGeo.type());
      bool hostEnIsInMacroCell = refElement.checkInside(macroGeo.local(baryCenter));
      if (hostEnIsInMacroCell) {
        lastIterator = macroCellIterator;
        auto macroIndex = coarseGridLeafIndexSet_.index(*macroCellIterator);
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

SubGridList::IdType SubGridList::getEnclosingMacroCellId(const HostEntityPointerType& hostEntityPointer) {
  // first check, whether we looked for this host entity already
  IdType hostEntityId = hostSpace_.gridPart().grid().globalIdSet().id(*hostEntityPointer);
  auto itFound = fineToCoarseMapID_.find(hostEntityId);
  if (itFound != fineToCoarseMapID_.end()) {
    // if so, return the index that was found last time
    return itFound->second;
  }
  static auto lastIterator = coarseSpace_.gridPart().grid().leafbegin<0>();
  const auto baryCenter = hostEntityPointer->geometry().center();
  auto macroCellIterator = lastIterator;
  for (; macroCellIterator != coarseSpace_.gridPart().grid().leafend<0>(); ++macroCellIterator) {
    const auto& macroGeo = macroCellIterator->geometry();
    const auto& refElement = CoarseRefElementType::general(macroGeo.type());

    bool hostEnIsInMacroCell = refElement.checkInside(macroGeo.local(baryCenter));
    if (hostEnIsInMacroCell) {
      lastIterator = macroCellIterator;
      IdType macroId = coarseSpace_.gridPart().grid().globalIdSet().id(*macroCellIterator);
      fineToCoarseMapID_[hostEntityId] = macroId;
      return macroId;
    }
  }
  // if we came this far, we did not find the matching enclosing coarse cell for the given
  // fine cell in [lastIterator, coarse grid end]. Start search from beginning
  for (macroCellIterator = coarseSpace_.gridPart().grid().leafbegin<0>(); macroCellIterator != lastIterator;
       ++macroCellIterator) {
    const auto& macroGeo = macroCellIterator->geometry();
    const auto& refElement = CoarseRefElementType::general(macroGeo.type());
    bool hostEnIsInMacroCell = refElement.checkInside(macroGeo.local(baryCenter));
    if (hostEnIsInMacroCell) {
      lastIterator = macroCellIterator;
      IdType macroId = coarseSpace_.gridPart().grid().globalIdSet().id(*macroCellIterator);
      fineToCoarseMapID_[hostEntityId] = macroId;
      return macroId;
    }
  }
  // if we came this far, we did not find an enclosing coarse cell at all, issue a warning
  // and return with error code
  //  DSC_LOG_DEBUG << "Warning: Host grid entity was not in any coarse grid cell!\n";
  assert(false);
  return -1;
}

void SubGridList::identifySubGrids() {
  DSC_PROFILER.startTiming("msfem.subgrid_list.identify");
  DSC_LOG_INFO << "Starting creation of subgrids." << std::endl << std::endl;

  // the fine grid part
  const HostGridPartType& hostGridPart = hostSpace_.gridPart();


  // -------- identify the entities that share a certain node -------

  // determine the entities that share a common global node with a given index
  // we need to iterate over the whole grid, not only from hostSpace_.begin() to
  // hostSpace_.end() for parallel runs!
  for (auto& hostEntity : DSC::viewRange(hostSpace_.gridPart().grid().leafView())) {
    int number_of_nodes_in_entity = hostEntity.count<HostGridType::dimension>();
    for (int i = 0; i < number_of_nodes_in_entity; ++i) {
      const HostNodePointer node = hostEntity.subEntity<HostGridType::dimension>(i);
      const auto global_index_node = hostGridPart.indexSet().index(*node);

      // make sure we don't access non-existing elements
      assert(entities_sharing_same_node_.size() > global_index_node);
      entities_sharing_same_node_[global_index_node].emplace_back(hostEntity);
    }
  }

  // the number of coarse grid entities (of codim 0).
  const auto number_of_coarse_grid_entities = specifier_.getNumOfCoarseEntities();
  const auto oversampling_strategy = specifier_.getOversamplingStrategy();

  if ((oversampling_strategy == 2) || (oversampling_strategy == 3)) {
    coarse_node_store_ = CoarseGridNodeStorageType(number_of_coarse_grid_entities, CoarseNodeVectorType());
    extended_coarse_node_store_ = CoarseGridNodeStorageType(number_of_coarse_grid_entities, CoarseNodeVectorType());
  }

  // ! ----------- create subgrids --------------------

  std::fill(enriched_.data(), enriched_.data() + enriched_.num_elements(), false);
  // loop to initialize subgrids (and to initialize the coarse node vector):
  // -----------------------------------------------------------

  // the fine grid (subgrid needs non-const ref)
  HostGridType& hostGrid = hostSpace_.gridPart().grid();

  //! creating SubGirds in parallel sections fail
  for (const auto& coarse_entity : coarseSpace_) {
    const auto coarse_index = coarseGridLeafIndexSet_.index(coarse_entity);
    // make sure we did not create a subgrid for the current coarse entity so far
    assert(subGridList_.find(coarse_index) == subGridList_.end());
    subGridList_[coarse_index] = std::make_shared<SubGridType>(hostGrid);
  }

  #ifdef _OPENMP
  #pragma omp parallel
  #endif
  {
  for (const auto& coarse_entity : threadIterators_) {
    // make sure we only create subgrids for interior coarse elements, not
    // for overlap or ghost elements
    assert(coarse_entity.partitionType() == Dune::InteriorEntity);
    const auto coarse_index = coarseGridLeafIndexSet_.index(coarse_entity);

    subgrid_id_to_base_coarse_entity_.insert(std::make_pair(coarse_index, std::move(coarse_entity.seed())));

    assert(subGridList_[coarse_index]);
    subGridList_[coarse_index]->createBegin();

    if ((oversampling_strategy == 2) || (oversampling_strategy == 3)) {
      assert(coarse_index >= 0 && coarse_index < coarse_node_store_.size() &&
             "Index set is not suitable for the current implementation!");
      for (int c = 0; c < coarse_entity.geometry().corners(); ++c) {
        coarse_node_store_[coarse_index].emplace_back(coarse_entity.geometry().corner(c));
        extended_coarse_node_store_[coarse_index].emplace_back(coarse_entity.geometry().corner(c));
      }
    }
  } // for
  } // omp region

  //  // -----------------------------------------------------------

  DSC_PROFILER.stopTiming("msfem.subgrid_list.identify",
                          DSC_CONFIG_GET("global.output_walltime", false));

  return;
}

void SubGridList::createSubGrids() {
  DSC_PROFILER.startTiming("msfem.subgrid_list.create");

  const HostGridLeafIndexSet& hostGridLeafIndexSet = hostSpace_.gridPart().grid().leafIndexSet();
  Fem::DomainDecomposedIteratorStorage< CommonTraits::GridPartType > hostSpaceIterators(hostSpace_.gridPart());
  hostSpaceIterators.update();

  // loop over all host entities and assign them to a macro cell
  #ifdef _OPENMP
  #pragma omp parallel
  #endif
  {
  for (const auto& host_entity : hostSpaceIterators) {
    // get the coarse-grid-father of host_entity (which is a maxlevel entity)...
    //    const HostEntityPointerType level_father_entity = DSG::make_father(coarseGridLeafIndexSet_,
    //            HostEntityPointerType(host_entity),
    //            specifier_.getLevelDifference());
    //// ... and its index

    const auto macroCellIndex = getEnclosingMacroCellIndex(host_entity);
    // add host_entity to the subgrid with the index 'macroCellIndex'
    subGridList_[macroCellIndex]->insertPartial(host_entity);
    // add the id of the subgrid to the vecor at position 'index of host grid element'
    fine_id_to_subgrid_ids_[hostGridLeafIndexSet.index(host_entity)].push_back(macroCellIndex);

    // check the neighbor entities and look if they belong to the same father
    // if yes, continue
    // if not, enrichment with 'n(T)' layers
    bool all_neighbors_have_same_father = true;
    const HostIntersectionIterator iend = hostGridPart_.iend(host_entity);
    for (HostIntersectionIterator iit = hostGridPart_.ibegin(host_entity);
         (iit != iend) && all_neighbors_have_same_father; ++iit) {
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
      auto layers = specifier_.getNoOfLayers(macroCellIndex);
      if (layers > 0) {
        DSC::Profiler::ScopedTiming enrichment_st("msfem.subgrid_list.enrichment");
        const HostEntityPointerType hep(host_entity);
        enrichment(hep, macroCellIndex, subGridList_[macroCellIndex], layers);
      }
    }

  } // for
  } // omp region

  DSC_PROFILER.stopTiming("msfem.subgrid_list.create",
                          DSC_CONFIG_GET("global.output_walltime", false));

  return;
}

template < class A, class B, class C>
void finalizeSubGrid(A& it, const B& coarseSpace_, const C& coarseGridLeafIndexSet_)
{
  auto& subGridIt = *it;
  subGridIt.second->createEnd();

  // report some infos about the subgrids if desired
  if (!false) {
    DSC_LOG_INFO << "Subgrid " << subGridIt.first << ":" << std::endl;
    subGridIt.second->report();
  }

  // error handling
  if (subGridIt.second->size(2) == 0) {
    DSC_LOG_ERROR << "Error." << std::endl << "Error. Created Subgrid with 0 nodes." << std::endl;

    for (const auto& coarse_it : coarseSpace_) {
      const auto index = coarseGridLeafIndexSet_.index(coarse_it);
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

void SubGridList::finalizeSubGrids() {
  DSC_PROFILER.startTiming("msfem.subgrid_list.create.finalize");

//#ifndef GNU
//#error bajh
//#else

//  auto func = std::bind(finalizeSubGrid, std::placeholders::_1, coarseSpace_, coarseGridLeafIndexSet_);
  auto func = [](SubGridStorageType::value_type& it){it.second->createEnd();};
//  func(*subGridList_.begin());
  __gnu_parallel::for_each(subGridList_.begin(),  subGridList_.end(), func);

  DSC_PROFILER.stopTiming("msfem.subgrid_list.create.finalize",
                          DSC_CONFIG_GET("global.output_walltime", false));

  return;
}

} // namespace MsFEM {
} // namespace Multiscale {
} // namespace Dune {
