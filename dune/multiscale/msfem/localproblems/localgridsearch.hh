#ifndef LOCALGRIDSEARCH_HH
#define LOCALGRIDSEARCH_HH

#include <dune/multiscale/msfem/msfem_traits.hh>
#include <dune/multiscale/msfem/localproblems/subgrid-list.hh>
#include <dune/stuff/grid/search.hh>
#include <dune/stuff/common/memory.hh>

namespace Dune {
namespace Multiscale {
namespace MsFEM {

template <class GridViewImp>
class LocalGridSearch
    : public DSG::EntitySearchBase<GridViewImp>
{
  typedef GridViewImp LocalGridViewType;
  typedef DSG::EntitySearchBase<LocalGridViewType> BaseType;
  typedef GenericReferenceElements< typename  BaseType::LocalCoordinateType::value_type,
                                             BaseType::LocalCoordinateType::dimension > RefElementType;
  typedef typename LocalGridViewType::template Codim< 0 >::Iterator IteratorType;

  typedef typename BaseType::EntityType::EntityPointer EntityPointer;

  typedef DSG::EntityInlevelSearch<LocalGridViewType> PerGridSearchType;
  typedef CommonTraits::DiscreteFunctionSpaceType CoarseGridSpaceType;
  typedef typename CoarseGridSpaceType::GridType::Traits::LeafIndexSet::IndexType IndexType;

public:
  typedef typename BaseType::EntityPointerVectorType EntityPointerVectorType;

  LocalGridSearch(const CoarseGridSpaceType& space, const LocalGridList& gridlist);

  template < class PointContainerType >
  EntityPointerVectorType operator() (const PointContainerType& points);

  const EntityPointer& current_coarse_pointer() const;
private:
  const CoarseGridSpaceType& coarse_space_;
  const LocalGridList& gridlist_;
  std::map<IndexType, std::unique_ptr<PerGridSearchType>> coarse_searches_;
};

template <class GridViewImp>
template <class PointContainerType>
typename LocalGridSearch<GridViewImp>::EntityPointerVectorType LocalGridSearch<GridViewImp>::operator()(const PointContainerType& points)
{
  auto count_nulls = [&](typename EntityPointerVectorType::value_type& ptr){return ptr==nullptr;};
  for(const auto& coarse_entity : coarse_space_) {
    const auto& localgrid = gridlist_.getSubGrid(coarse_entity);
    const auto& index_set = coarse_space_.gridPart().grid().leafIndexSet();
    const auto index = index_set.index(coarse_entity);
    auto& current_search_ptr = coarse_searches_[index];
    if(!current_search_ptr)
      current_search_ptr = DSC::make_unique<PerGridSearchType>(localgrid.leafView());
    auto entity_ptrs = current_search_ptr->operator()(points);
    auto non_nulls = std::count_if(entity_ptrs.begin(), entity_ptrs.end(), count_nulls);
    if(non_nulls == 0)
      return entity_ptrs;
    if(non_nulls < long(entity_ptrs.size()))
      DUNE_THROW(InvalidStateException, "partial coverage not supported atm");
  }
  DUNE_THROW(InvalidStateException, "local grid search failed");
}

template <class GridViewImp>
LocalGridSearch<GridViewImp>::LocalGridSearch(const CoarseGridSpaceType& coarse_space, const LocalGridList& gridlist)
  : coarse_space_(coarse_space)
  , gridlist_(gridlist)
{}

} // namespace MsFEM {
} // namespace Multiscale {
} // namespace Dune {

#endif // LOCALGRIDSEARCH_HH
