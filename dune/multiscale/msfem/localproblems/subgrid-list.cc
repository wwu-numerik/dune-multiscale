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
#include <dune/common/parallel/mpihelper.hh>
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
  , coarseGridLeafIndexSet_(coarseSpace_.grid_view()->grid().leafIndexSet()) {
  DSC::Profiler::ScopedTiming st("msfem.subgrid_list.ctor");
  BOOST_ASSERT_MSG(DSC_CONFIG.hasSub("grids"), "Parameter tree needs to have 'grids' subtree!");
  const int dim_world = LocalGridType::dimensionworld;

  const DSC::ExtendedParameterTree gridParameterTree(DSC_CONFIG.sub("grids"));
  std::vector<int> micro_per_macro = gridParameterTree.getVector("micro_cells_per_macrocell_dim", 8, dim_world);
  const auto oversampling_layer = DSC_CONFIG_GET("msfem.oversampling_layers", 0);

  typedef StructuredGridFactory<LocalGridType> FactoryType;
  const auto coarse_dimensions = DSG::dimensions<CommonTraits::GridType>(coarseSpace_.grid_view()->grid());

  for (const auto& coarse_entity : DSC::viewRange(*coarseSpace_.grid_view())) {
    // make sure we only create subgrids for interior coarse elements, not
    // for overlap or ghost elements
    assert(coarse_entity.partitionType() == Dune::InteriorEntity);
    const auto coarse_index = coarseGridLeafIndexSet_.index(coarse_entity);
    // make sure we did not create a subgrid for the current coarse entity so far
    assert(subGridList_.find(coarse_index) == subGridList_.end());

    const auto dimensions = DSG::dimensions<CommonTraits::GridType>(coarse_entity);
    typedef FieldVector<typename LocalGridType::ctype, dim_world> CoordType;
    CoordType lowerLeft(0);
    CoordType upperRight(0);
    array<unsigned int, dim_world> elemens;
    for (const auto i : DSC::valueRange(dim_world)) {
      const auto min = dimensions.coord_limits[i].min();
      const auto max = dimensions.coord_limits[i].max();
      auto localMin = coarse_dimensions.coord_limits[i].min();
      auto localMax = coarse_dimensions.coord_limits[i].max();
      auto comm = MPIHelper::getCollectiveCommunication();
      const auto coarse_min = comm.min(localMin);
      const auto coarse_max = comm.max(localMax);
      const auto delta = (max - min) / double(micro_per_macro[i]);
      lowerLeft[i] = std::max(min - (oversampling_layer * delta), coarse_min);
      upperRight[i] = std::min(max + (oversampling_layer * delta), coarse_max);
      int smaller = ((min - (oversampling_layer * delta)) < coarse_min);
      int bigger = ((max + (oversampling_layer * delta)) > coarse_max);
      elemens[i] = micro_per_macro[i] + ((!smaller + !bigger) * oversampling_layer);
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
