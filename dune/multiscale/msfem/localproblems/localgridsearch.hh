#ifndef DUNE_MULTISCALE_MSFEM_LOCALGRIDSEARCH_HH
#define DUNE_MULTISCALE_MSFEM_LOCALGRIDSEARCH_HH

#include <dune/multiscale/msfem/msfem_traits.hh>
#include <dune/stuff/grid/search.hh>


namespace Dune {
namespace Multiscale {
namespace MsFEM {

class LocalGridList;

//! given a Localgridlist, facilitate searching for evaluation points in a pseudo-hierachical manner
class LocalGridSearch : public DSG::EntitySearchBase<MsFEMTraits::LocalGridViewType> {
  typedef MsFEMTraits::LocalGridViewType LocalGridViewType;
  typedef DSG::EntitySearchBase<LocalGridViewType> BaseType;
  typedef typename LocalGridViewType::template Codim<0>::Iterator IteratorType;

  typedef typename BaseType::EntityType::EntityPointer EntityPointer;

  typedef DSG::EntityInlevelSearch<LocalGridViewType> PerGridSearchType;
  typedef CommonTraits::DiscreteFunctionSpaceType CoarseGridSpaceType;
  typedef typename CoarseGridSpaceType::GridViewType::Grid::Traits::LeafIndexSet::IndexType IndexType;
  typedef typename CoarseGridSpaceType::EntityType::EntityPointer CoarseEntityPointerType;
  typedef std::vector<CoarseGridSpaceType::DomainType> PointContainerType;
  typedef PointContainerType::const_iterator PointIterator;

public:
  typedef typename BaseType::EntityPointerVectorType EntityPointerVectorType;

  LocalGridSearch(const CoarseGridSpaceType& space, const LocalGridList& gridlist);

  EntityPointerVectorType operator()(const PointContainerType& points);

  const CoarseEntityPointerType& current_coarse_pointer() const;

  bool covers_strict(const typename CoarseGridSpaceType::EntityType& coarse_entity, const PointIterator first,
                     const PointIterator last);

private:
  const CoarseGridSpaceType& coarse_space_;
  const LocalGridList& gridlist_;
  std::map<IndexType, std::unique_ptr<PerGridSearchType>> coarse_searches_;
  std::unique_ptr<CoarseEntityPointerType> current_coarse_pointer_;
};


} // namespace MsFEM {
} // namespace Multiscale {
} // namespace Dune {

#endif // DUNE_MULTISCALE_MSFEM_LOCALGRIDSEARCH_HH
