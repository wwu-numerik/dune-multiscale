#ifndef DUNE_MULTISCALE_MSFEM_LOCALGRIDSEARCH_HH
#define DUNE_MULTISCALE_MSFEM_LOCALGRIDSEARCH_HH

#include <dune/multiscale/msfem/msfem_traits.hh>
#include <dune/stuff/grid/search.hh>

namespace Dune {
namespace Multiscale {

class LocalGridList;

//! given a Localgridlist, facilitate searching for evaluation points in a pseudo-hierachical manner
class LocalGridSearch : public DSG::EntitySearchBase<MsFEMTraits::LocalGridViewType> {
  typedef MsFEMTraits::LocalGridViewType LocalGridViewType;
  typedef DSG::EntitySearchBase<LocalGridViewType> BaseType;

  typedef DSG::EntityInlevelSearch<LocalGridViewType> PerGridSearchType;
  typedef CommonTraits::SpaceType CoarseGridSpaceType;
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
  CommonTraits::InteriorGridViewType static_view_;
  typedef typename CommonTraits::InteriorGridViewType::template Codim<0>::Iterator InteriorIteratorType;
  std::unique_ptr<InteriorIteratorType> static_iterator_;
};

} // namespace Multiscale {
} // namespace Dune {

#endif // DUNE_MULTISCALE_MSFEM_LOCALGRIDSEARCH_HH
