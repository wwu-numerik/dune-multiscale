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


#include "subgrid-list.hh"

namespace Dune {
namespace Multiscale {
namespace MsFEM {

LocalGridList::LocalGridList(const CommonTraits::DiscreteFunctionSpaceType& coarseSpace)
  : coarseSpace_(coarseSpace)
  , coarseGridLeafIndexSet_(coarseSpace_.gridPart().grid().leafIndexSet()) {
  DSC::Profiler::ScopedTiming st("msfem.subgrid_list.ctor");
  const auto micro_per_macro = DSC_CONFIG_GET("msfem.micro_cells_per_macrocell_dim", 8);
  const auto oversampling_layer = DSC_CONFIG_GET("msfem.oversampling_layers", 0);

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
    array<unsigned int, dim_world> elemens;
    for (const auto i : DSC::valueRange(dim_world)) {
      elemens[i] = micro_per_macro + (2 * oversampling_layer);
      const auto min = dimensions.coord_limits[i].min();
      const auto max = dimensions.coord_limits[i].max();
      const auto coarse_min = coarse_dimensions.coord_limits[i].min();
      const auto coarse_max = coarse_dimensions.coord_limits[i].max();
      const auto delta = (max - min) / double(micro_per_macro);
      lowerLeft[i] = std::max(min - (oversampling_layer * delta), coarse_min);
      upperRight[i] = std::min(max + (oversampling_layer * delta), coarse_max);
    }
    subGridList_[coarse_index] = FactoryType::createCubeGrid(lowerLeft, upperRight, elemens);
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

bool LocalGridList::covers_strict(const CoarseEntityType& coarse_entity,
                                  const MsFEMTraits::LocalEntityType& local_entity) {
  const auto& reference_element = Stuff::Grid::reference_element(coarse_entity);
  const auto& coarse_geometry = coarse_entity.geometry();
  for (const auto i : DSC::valueRange(local_entity.geometry().corners())) {
    if (!reference_element.checkInside(coarse_geometry.local(local_entity.geometry().corner(i))))
      return false;
  }
  return true;
}

bool LocalGridList::covers(const CoarseEntityType& coarse_entity, const MsFEMTraits::LocalEntityType& local_entity) {
  const auto& center = local_entity.geometry().center();
  const auto& reference_element = Stuff::Grid::reference_element(coarse_entity);
  return reference_element.checkInside(coarse_entity.geometry().local(center));
}

} // namespace MsFEM {
} // namespace Multiscale {
} // namespace Dune {
