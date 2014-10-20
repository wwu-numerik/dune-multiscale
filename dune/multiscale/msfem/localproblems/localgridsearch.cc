#include <config.h>

#include "localgridsearch.hh"

#include <dune/multiscale/msfem/localproblems/localgridlist.hh>
#include <dune/stuff/common/memory.hh>
#include <dune/stuff/common/algorithm.hh>
#include <dune/grid/common/gridenums.hh>

Dune::Multiscale::LocalGridSearch::EntityPointerVectorType Dune::Multiscale::LocalGridSearch::
operator()(const PointContainerType& points) {
  typedef typename EntityPointerVectorType::value_type EPV;
  const auto is_null = [&](const EPV& ptr) { return ptr == nullptr; };
  const auto not_null = [&](const EPV& ptr) { return ptr != nullptr; };

  // only iterate over inner (non-overlap) entities
  static const auto view =
      coarse_space_.grid_view().grid().template leafGridView<PartitionIteratorType::InteriorBorder_Partition>();

  static auto it = view.template begin<0>();
  const auto end = view.template end<0>();
  EntityPointerVectorType ret_entities(points.size());
  int steps = 0;
  bool did_cover = false;
  auto null_count = points.size();
  while (null_count) {
    assert(it != end);
    const auto& coarse_entity = *it;
    const auto& localgrid = gridlist_.getSubGrid(coarse_entity);
    const auto& index_set = view.grid().leafIndexSet();
    const auto index = index_set.index(coarse_entity);
    current_coarse_pointer_ = DSC::make_unique<CoarseEntityPointerType>(coarse_entity);
    auto& current_search_ptr = coarse_searches_[index];
    if (current_search_ptr == nullptr)
      current_search_ptr = DSC::make_unique<PerGridSearchType>(localgrid.leafGridView());
    did_cover = covers_strict(coarse_entity, points.begin(), points.end());
    if (did_cover) {
      auto first_null = std::find(ret_entities.begin(), ret_entities.end(), nullptr);
      auto entity_ptrs = current_search_ptr->operator()(points);
      DSC::move_if(entity_ptrs.begin(), entity_ptrs.end(), first_null, not_null);
      null_count = std::count_if(ret_entities.begin(), ret_entities.end(), is_null);
    }
    if (++it == end) {
      if (++steps < view.size(0))
        it = view.template begin<0>();
      else {
        DUNE_THROW(InvalidStateException, "local grid search failed ");
      }
    }
  }
  return ret_entities;
}

bool Dune::Multiscale::LocalGridSearch::covers_strict(const CoarseGridSpaceType::EntityType& coarse_entity,
                                                      const Dune::Multiscale::LocalGridSearch::PointIterator first,
                                                      const Dune::Multiscale::LocalGridSearch::PointIterator last) {
  const auto& reference_element = Stuff::Grid::reference_element(coarse_entity);
  const auto& coarse_geometry = coarse_entity.geometry();
  for (auto it = first; it != last; ++it) {
    if (!reference_element.checkInside(coarse_geometry.local(*it)))
      return false;
  }
  return true;
}

Dune::Multiscale::LocalGridSearch::LocalGridSearch(const CoarseGridSpaceType& coarse_space,
                                                   const Dune::Multiscale::LocalGridList& gridlist)
  : coarse_space_(coarse_space)
  , gridlist_(gridlist) {}

const Dune::Multiscale::LocalGridSearch::CoarseEntityPointerType&
Dune::Multiscale::LocalGridSearch::current_coarse_pointer() const {
  assert(current_coarse_pointer_);
  return *current_coarse_pointer_;
}
