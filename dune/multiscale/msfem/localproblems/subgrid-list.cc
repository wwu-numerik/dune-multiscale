#include <config.h>
#include <assert.h>
#include <boost/assert.hpp>
#include <boost/multi_array/multi_array_ref.hpp>
#include <dune/common/exceptions.hh>
#include <dune/fem/misc/mpimanager.hh>
#include <dune/stuff/common/logging.hh>
#include <dune/stuff/common/profiler.hh>
#include <dune/stuff/common/ranges.hh>
#include <dune/stuff/common/float_cmp.hh>
#include <dune/stuff/grid/structuredgridfactory.hh>
#include <dune/stuff/grid/information.hh>
#include <algorithm>
#include <iterator>
#include <ostream>
#include <utility>
#include <memory>
#include <dune/multiscale/tools/misc.hh>

#include <Eigen/Core>

#include "subgrid-list.hh"

namespace Dune {
namespace Multiscale {
namespace MsFEM {

LocalGridList::LocalGridList(MsFEMTraits::MacroMicroGridSpecifierType& specifier, bool silent /*= true*/)
  : coarseSpace_(specifier.coarseSpace())
  , specifier_(specifier)
  , silent_(silent)
  , coarseGridLeafIndexSet_(coarseSpace_.gridPart().grid().leafIndexSet())
  #if 0 // LOD only
  , fineToCoarseMap_(Fem::MPIManager::size())
  #endif //0 // LOD only
{
  DSC::Profiler::ScopedTiming st("msfem.subgrid_list.ctor");
  const auto oversampling_strategy = DSC_CONFIG_GET("msfem.oversampling_strategy", 1);
  const auto micro_per_macro = DSC_CONFIG_GET("msfem.micro_cells_per_macrocell_dim", 8);
  const auto oversampling_layer = DSC_CONFIG_GET("msfem.oversampling_layers", 0);
#if 0 // LOD only
  const auto number_of_coarse_grid_entities = specifier_.getNumOfCoarseEntities();
  if ((oversampling_strategy == 2) || (oversampling_strategy == 3)) {
    coarse_node_store_ = CoarseGridNodeStorageType(number_of_coarse_grid_entities, CoarseNodeVectorType());
    extended_coarse_node_store_ = CoarseGridNodeStorageType(number_of_coarse_grid_entities, CoarseNodeVectorType());
  }
#endif // 0 // LOD only
  typedef StructuredGridFactory<LocalGridType> FactoryType;
  const auto coarse_dimensions = DSG::dimensions<CommonTraits::GridType>(coarseSpace_.gridPart().grid());

  for (const auto& coarse_entity : coarseSpace_) {
    // make sure we only create subgrids for interior coarse elements, not
    // for overlap or ghost elements
    assert(coarse_entity.partitionType() == Dune::InteriorEntity);
    const auto coarse_index = coarseGridLeafIndexSet_.index(coarse_entity);
    // make sure we did not create a subgrid for the current coarse entity so far
    assert(subGridList_.find(coarse_index) == subGridList_.end());

    const auto dimensions = DSG::dimensions<CommonTraits::GridType>(coarse_entity);
    const int dim_world = LocalGridType::dimensionworld;
    typedef FieldVector<typename LocalGridType::ctype, dim_world> CoordType;
    CoordType lowerLeft(0);
    CoordType upperRight(0);
    array<unsigned int,dim_world> elemens;
    for(const auto i : DSC::valueRange(dim_world))
    {
      elemens[i] = micro_per_macro + (2 * oversampling_layer);
      const auto min = dimensions.coord_limits[i].min();
      const auto max = dimensions.coord_limits[i].max();
      const auto coarse_min = coarse_dimensions.coord_limits[i].min();
      const auto coarse_max = coarse_dimensions.coord_limits[i].max();
      const auto delta = (max - min) / double(micro_per_macro);
      lowerLeft[i] = std::max(min - (oversampling_layer * delta), coarse_min);
      upperRight[i] = std::min(max + (oversampling_layer * delta), coarse_max);
    }

    boost::format sp("LocalGrid %d from (%f,%f) to (%f,%f) created.\n");
    DSC_LOG_INFO << sp % coarse_index % lowerLeft[0] % lowerLeft[1] % upperRight[0] % upperRight[1];
    subGridList_[coarse_index] = FactoryType::createCubeGrid(lowerLeft, upperRight, elemens);

    if ((oversampling_strategy == 2) || (oversampling_strategy == 3)) {
      assert(coarse_index >= 0 && coarse_index < coarse_node_store_.size() &&
             "Index set is not suitable for the current implementation!");
      for (int c = 0; c < coarse_entity.geometry().corners(); ++c) {
        coarse_node_store_[coarse_index].emplace_back(coarse_entity.geometry().corner(c));
#if 0 // LOD only
        extended_coarse_node_store_[coarse_index].emplace_back(coarse_entity.geometry().corner(c));
#endif // 0 // LOD only
      }
    }
#if 0 // LOD only
    subgrid_id_to_base_coarse_entity_.insert(std::make_pair(coarse_index, std::move(coarse_entity.seed())));
#endif // 0 // LOD only
  }
}

MsFEMTraits::LocalGridType& LocalGridList::getSubGrid(IndexType coarseCellIndex) {
  auto found = subGridList_.find(coarseCellIndex);
  BOOST_ASSERT_MSG(found != subGridList_.end(), "There is no subgrid for the index you provided!");
  assert(found->second);
  return *(found->second);
} // getSubGrid


const MsFEMTraits::LocalGridType& LocalGridList::getSubGrid(IndexType coarseCellIndex) const {
  auto found = subGridList_.find(coarseCellIndex);
  BOOST_ASSERT_MSG(found != subGridList_.end(), "There is no subgrid for the index you provided!");
  assert(found->second);
  return *(found->second);
} // getSubGrid

const MsFEMTraits::LocalGridType& LocalGridList::getSubGrid(const CoarseEntityType& entity) const {
  const auto index = coarseGridLeafIndexSet_.index(entity);
  return getSubGrid(index);
} // getSubGrid

MsFEMTraits::LocalGridType& LocalGridList::getSubGrid(const CoarseEntityType& entity) {
  const auto index = coarseGridLeafIndexSet_.index(entity);
  return getSubGrid(index);
} // getSubGrid

// get number of sub grids
std::size_t LocalGridList::size() const { return subGridList_.size(); }

MsFEM::MsFEMTraits::LocalGridPartType LocalGridList::gridPart(IndexType i) { return LocalGridPartType(getSubGrid(i)); }

bool LocalGridList::covers(const CoarseEntityType &coarse_entity, const LocalEntityType &local_entity) {
  const auto& center = local_entity.geometry().center();
  const auto& coarse_geo = coarse_entity.geometry();
  const auto center_local = coarse_geo.local(center);
  const auto& reference_element = Stuff::Grid::reference_element(coarse_entity);
  return reference_element.checkInside(center_local);
}

const LocalGridList::CoarseNodeVectorType& LocalGridList::getCoarseNodeVector(IndexType i) const {
  if (DSC_CONFIG_GET("msfem.oversampling_strategy", 1) == 1)
    DUNE_THROW(Dune::InvalidStateException, "Method 'getCoarseNodeVector' of class 'LocalGridList' should not be used in\
                combination with oversampling strategy 1. Check your implementation!");

  if (i >= specifier_.getNumOfCoarseEntities()) {
    DUNE_THROW(Dune::RangeError, "Error. LocalGrid-Index too large.");
  }
  return coarse_node_store_[i];
} // getSubGrid

#if 0 // LOD only
// given the id of a subgrid, return the entity seed for the 'base coarse entity'
// (i.e. the coarse entity that the subgrid was constructed from by enrichment )
const LocalGridList::CoarseEntitySeedType &LocalGridList::get_coarse_entity_seed(std::size_t i) const {
  // the following returns the mapped element for index i if present,
  // if not, an out-of-range exception is thrown
  assert(false);//need to eliminate narrowing conversion
  return subgrid_id_to_base_coarse_entity_.at(i);
}

const LocalGridList::CoarseNodeVectorType& LocalGridList::getExtendedCoarseNodeVector(IndexType i) const {
  if ((DSC_CONFIG_GET("msfem.oversampling_strategy", 1) == 1) || (DSC_CONFIG_GET("msfem.oversampling_strategy", 1) == 2))
    DUNE_THROW(Dune::InvalidStateException,
               "Method 'getExtendendCoarseNodeVector' of class 'LocalGridList' should not be used in\
                combination with oversampling strategy 1 or 2. Check your implementation!");

  if (i >= specifier_.getNumOfCoarseEntities()) {
    DUNE_THROW(Dune::RangeError, "Error. LocalGrid-Index too large.");
  }
  return extended_coarse_node_store_[i];
} // getSubGrid
#endif // 0 // LOD only

} // namespace MsFEM {
} // namespace Multiscale {
} // namespace Dune {
