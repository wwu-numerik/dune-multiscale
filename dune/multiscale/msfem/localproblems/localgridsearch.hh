#ifndef LOCALGRIDSEARCH_HH
#define LOCALGRIDSEARCH_HH

#include <dune/multiscale/msfem/msfem_traits.hh>
#include <dune/multiscale/msfem/localproblems/subgrid-list.hh>
#include <dune/stuff/grid/search.hh>

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

public:
  typedef typename BaseType::EntityPointerVectorType EntityPointerVectorType;

  LocalGridSearch(const LocalGridList& gridlist);

  template < class PointContainerType >
  EntityPointerVectorType operator() (const PointContainerType& points);

  const EntityPointer& current_coarse_pointer() const;
private:
  const LocalGridList& gridlist_;
};

template <class GridViewImp>
template <class PointContainerType>
typename LocalGridSearch<GridViewImp>::EntityPointerVectorType LocalGridSearch<GridViewImp>::operator()(const PointContainerType& points)
{

}

template <class GridViewImp>
LocalGridSearch<GridViewImp>::LocalGridSearch(const LocalGridList& gridlist)
  : gridlist_(gridlist)
{}

} // namespace MsFEM {
} // namespace Multiscale {
} // namespace Dune {

#endif // LOCALGRIDSEARCH_HH
