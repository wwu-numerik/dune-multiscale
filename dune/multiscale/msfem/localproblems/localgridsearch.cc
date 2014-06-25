#include <config.h>

#include "localgridsearch.hh"

#include <dune/multiscale/msfem/localproblems/subgrid-list.hh>
#include <dune/stuff/common/memory.hh>
#include <dune/grid/common/gridenums.hh>

Dune::Multiscale::MsFEM::LocalGridSearch::EntityPointerVectorType Dune::Multiscale::MsFEM::LocalGridSearch::operator()(const PointContainerType &points) {
  const auto count_nulls = [&](const typename EntityPointerVectorType::value_type& ptr) { return ptr == nullptr; };
  //! \TODO potential speedup by caching last coarse_entity position instead fo restarting at front
  // only iterate over inner (non-overlap) entities
  const auto view = coarse_space_.grid_view()->grid().template leafGridView<PartitionIteratorType::Interior_Partition>();

  static auto it = view.template begin< 0 >();
  const auto end = view.template end< 0 >();
  int steps = 0;
  while(true) {
    const auto& coarse_entity = *it;
    const auto& localgrid = gridlist_.getSubGrid(coarse_entity);
    const auto& index_set = view.grid().leafIndexSet();
    const auto index = index_set.index(coarse_entity);
    current_coarse_pointer_ = DSC::make_unique<CoarseEntityPointerType>(coarse_entity);
    auto& current_search_ptr = coarse_searches_[index];
    if (current_search_ptr == nullptr)
      current_search_ptr = DSC::make_unique<PerGridSearchType>(localgrid.leafGridView());
    if (covers_strict(coarse_entity, points.begin(), points.end())) {
      auto entity_ptrs = current_search_ptr->operator()(points);
      const auto null_count = std::count_if(entity_ptrs.begin(), entity_ptrs.end(), count_nulls);
      if (null_count == 0)
        return entity_ptrs;
    }
    if (++it == end) {
      if (++steps < view.size(0))
        it = view.template begin< 0 >();
      else
        DUNE_THROW(InvalidStateException, "local grid search failed");
    }
  }
}


bool Dune::Multiscale::MsFEM::LocalGridSearch::covers_strict(const CoarseGridSpaceType::EntityType &coarse_entity, const Dune::Multiscale::MsFEM::LocalGridSearch::PointIterator first, const Dune::Multiscale::MsFEM::LocalGridSearch::PointIterator last) {
  const auto& reference_element = Stuff::Grid::reference_element(coarse_entity);
  const auto& coarse_geometry = coarse_entity.geometry();
  for (auto it = first; it != last; ++it) {
    if (!reference_element.checkInside(coarse_geometry.local(*it)))
      return false;
  }
  return true;
}


Dune::Multiscale::MsFEM::LocalGridSearch::LocalGridSearch(const CoarseGridSpaceType &coarse_space, const Dune::Multiscale::MsFEM::LocalGridList &gridlist)
  : coarse_space_(coarse_space)
  , gridlist_(gridlist) {}


const Dune::Multiscale::MsFEM::LocalGridSearch::CoarseEntityPointerType &Dune::Multiscale::MsFEM::LocalGridSearch::current_coarse_pointer() const {
  assert(current_coarse_pointer_);
  return *current_coarse_pointer_;
}
