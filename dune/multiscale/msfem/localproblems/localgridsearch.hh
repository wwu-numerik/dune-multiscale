#ifndef DUNE_MULTISCALE_MSFEM_LOCALGRIDSEARCH_HH
#define DUNE_MULTISCALE_MSFEM_LOCALGRIDSEARCH_HH

#include <dune/multiscale/msfem/msfem_traits.hh>
#include <dune/stuff/grid/search.hh>

namespace Dune {
namespace Multiscale {

class LocalGridList;

//! given a Localgridlist, facilitate searching for evaluation points in a pseudo-hierachical manner
class LocalGridSearch : public DSG::EntitySearchBase<MsFEMTraits::LocalGridViewType>
{
  typedef DSG::EntitySearchBase<MsFEMTraits::LocalGridViewType> BaseType;
  typedef DSG::EntityInlevelSearch<MsFEMTraits::LocalGridViewType> PerGridSearchType;
  typedef typename CommonTraits::SpaceType::GridViewType::Grid::Traits::LeafIndexSet::IndexType IndexType;
  typedef typename CommonTraits::SpaceType::EntityType::EntityPointer CoarseEntityPointerType;
  typedef std::vector<CommonTraits::DomainType> PointContainerType;
  typedef PointContainerType::const_iterator PointIterator;

public:
  typedef typename BaseType::EntityVectorType EntityVectorType;

  LocalGridSearch(const CommonTraits::SpaceType& space, const LocalGridList& gridlist);
  LocalGridSearch(const LocalGridSearch& other);

  EntityVectorType operator()(const PointContainerType& points);

  const MsFEMTraits::CoarseEntityType& current_coarse_pointer() const;

  bool covers_strict(const CommonTraits::SpaceType::EntityType& coarse_entity,
                     const PointIterator first,
                     const PointIterator last);

private:
  const CommonTraits::SpaceType& coarse_space_;
  const LocalGridList& gridlist_;
  std::map<IndexType, std::unique_ptr<PerGridSearchType>> coarse_searches_;
  std::unique_ptr<MsFEMTraits::CoarseEntityType> current_coarse_entity_;
  CommonTraits::InteriorGridViewType static_view_;
  typedef typename CommonTraits::InteriorGridViewType::template Codim<0>::Iterator InteriorIteratorType;
  std::unique_ptr<InteriorIteratorType> static_iterator_;
};

} // namespace Multiscale {
} // namespace Dune {

#endif // DUNE_MULTISCALE_MSFEM_LOCALGRIDSEARCH_HH
