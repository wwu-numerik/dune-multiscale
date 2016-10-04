#include <config.h>

#include <algorithm>
#include <assert.h>
#include <boost/assert.hpp>
#include <boost/multi_array/multi_array_ref.hpp>
#include <dune/common/exceptions.hh>
#include <dune/multiscale/common/mygridfactory.hh>
#include <dune/multiscale/problems/selector.hh>
#include <dune/multiscale/tools/misc.hh>
#include <dune/stuff/grid/information.hh>
#include <dune/stuff/grid/structuredgridfactory.hh>
#include <dune/xt/common/float_cmp.hh>
#include <dune/xt/common/logging.hh>
#include <dune/xt/common/ranges.hh>
#include <dune/xt/common/timings.hh>
#include <iterator>
#include <memory>
#include <ostream>
#include <utility>

#include "localgridlist.hh"

namespace Dune {
namespace Multiscale {

LocalGridList::LocalGridList(const Problem::ProblemContainer& problem, const CommonTraits::SpaceType& coarseSpace)
  : coarseSpace_(coarseSpace)
  , coarseGridLeafIndexSet_(coarseSpace_.grid_view().grid().leafIndexSet()) {
  Dune::XT::Common::ScopedTiming algo("msfem.local_grids");
  BOOST_ASSERT_MSG(DXTC_CONFIG.has_sub("grids"), "Parameter tree needs to have 'grids' subtree!");
  constexpr auto dim_world = MsFEMTraits::LocalGridType::dimensionworld;

  const auto gridParameterTree = DXTC_CONFIG.sub("grids");
  const auto micro_per_macro = gridParameterTree.get<CommonTraits::DomainType>("micro_cells_per_macrocell_dim",
                                                                               CommonTraits::DomainType(8), dim_world);
  const auto oversampling_layer = problem.config().get("msfem.oversampling_layers", 0);

  typedef MyGridFactory<MsFEMTraits::LocalGridType> FactoryType;
  const auto& gridCorners = problem.getModelData().gridCorners();
  const auto globalLowerLeft = gridCorners.first;
  const auto globalUpperRight = gridCorners.second;

  const auto interior = interior_border_view(coarseSpace_);
  for (const auto& coarse_entity : elements(interior)) {
    // make sure we only create subgrids for interior coarse elements, not
    // for overlap or ghost elements
    assert(coarse_entity.partitionType() == Dune::InteriorEntity);
    const auto coarse_index = coarseGridLeafIndexSet_.index(coarse_entity);
    // make sure we did not create a subgrid for the current coarse entity so far
    assert(subGridList_.find(coarse_index) == subGridList_.end());

    const auto dimensions = DSG::dimensions<CommonTraits::GridType::LeafGridView>(coarse_entity);
    typedef FieldVector<typename MsFEMTraits::LocalGridType::ctype, dim_world> CoordType;
    CoordType lowerLeft(0);
    CoordType upperRight(0);
    array<unsigned int, dim_world> elements_per_dim;
    for (const auto i : Dune::XT::Common::value_range(dim_world)) {
      const auto min = dimensions.coord_limits[i].min();
      const auto max = dimensions.coord_limits[i].max();

      const auto delta = (max - min) / double(micro_per_macro[i]);
      lowerLeft[i] = std::max(min - (oversampling_layer * delta), globalLowerLeft[i]);
      upperRight[i] = std::min(max + (oversampling_layer * delta), globalUpperRight[i]);
      const bool exceeds_macro_grid_left = ((min - (oversampling_layer * delta)) < globalLowerLeft[i]);
      const bool exceeds_macro_grid_right = ((max + (oversampling_layer * delta)) > globalUpperRight[i]);
      elements_per_dim[i] = micro_per_macro[i] + ((int(!exceeds_macro_grid_left) + int(!exceeds_macro_grid_right)) * oversampling_layer);
    }
    subGridList_[coarse_index] = FactoryType::createLocalGrid(lowerLeft, upperRight, elements_per_dim);
  }
}

MsFEMTraits::LocalGridType& LocalGridList::getSubGrid(IndexType coarseCellIndex) {
  auto found = subGridList_.find(coarseCellIndex);
  BOOST_ASSERT_MSG(found != subGridList_.end(), "There is no subgrid for the index you provided!");
  assert(found->second);
  return *(found->second);
} // getSubGrid

const MsFEMTraits::LocalGridType& LocalGridList::getSubGrid(IndexType coarseCellIndex) const {
  const auto found = subGridList_.find(coarseCellIndex);
  BOOST_ASSERT_MSG(found != subGridList_.end(), "There is no subgrid for the index you provided!");
  assert(found->second);
  return *(found->second);
} // getSubGrid

const MsFEMTraits::LocalGridType& LocalGridList::getSubGrid(const MsFEMTraits::CoarseEntityType& entity) const {
  BOOST_ASSERT_MSG((entity.partitionType() == Dune::InteriorEntity), "Subgrids exist only for interior entities!");
  const auto index = coarseGridLeafIndexSet_.index(entity);
  return getSubGrid(index);
} // getSubGrid

MsFEMTraits::LocalGridType& LocalGridList::getSubGrid(const MsFEMTraits::CoarseEntityType& entity) {
  BOOST_ASSERT_MSG((entity.partitionType() == Dune::InteriorEntity), "Subgrids exist only for interior entities!");
  const auto index = coarseGridLeafIndexSet_.index(entity);
  return getSubGrid(index);
} // getSubGrid

// get number of sub grids
std::size_t LocalGridList::size() const { return subGridList_.size(); }

bool LocalGridList::covers_strict(const MsFEMTraits::CoarseEntityType& coarse_entity,
                                  const MsFEMTraits::LocalEntityType& local_entity) {
  return covers_strict(coarse_entity, local_entity.geometry());
}

bool LocalGridList::covers(const MsFEMTraits::CoarseEntityType& coarse_entity,
                           const MsFEMTraits::LocalEntityType& local_entity) {
  const auto& center = local_entity.geometry().center();
  const auto& reference_element = Stuff::Grid::reference_element(coarse_entity);
  return reference_element.checkInside(coarse_entity.geometry().local(center));
}

} // namespace Multiscale {
} // namespace Dune {
