// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef SUBGRIDLIST_HH
#define SUBGRIDLIST_HH

#include <dune/multiscale/common/traits.hh>
#include <dune/multiscale/msfem/msfem_traits.hh>

#include <boost/noncopyable.hpp>
#include <dune/common/deprecated.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/shared_ptr.hh>
#include <dune/stuff/grid/entity.hh>

#include <cstddef>
#include <map>
#include <memory>
#include <vector>

namespace Dune {
namespace Multiscale {


//! container for cell problem subgrids
class LocalGridList : public boost::noncopyable {
  typedef typename MsFEMTraits::LocalGridDiscreteFunctionType LocalGridDiscreteFunctionType;
  typedef typename MsFEMTraits::LocalSpaceType LocalSpaceType;
  typedef typename MsFEMTraits::LocalGridType LocalGridType;

  typedef typename CommonTraits::GridType::Traits::LeafIndexSet LeafIndexSet;
  typedef typename LeafIndexSet::IndexType IndexType;

  typedef typename MsFEMTraits::CoarseEntityType CoarseEntityType;
  typedef typename CoarseEntityType::EntitySeed CoarseEntitySeedType;
  typedef typename LeafIndexSet::IndexType CoarseEntityIndexType;
  typedef typename CommonTraits::DomainType DomainType;

public:
  typedef std::vector<DomainType> CoarseNodeVectorType;

  LocalGridList(const CommonTraits::SpaceType& coarseSpace);

  // private:
  LocalGridType& getSubGrid(IndexType i);
  const LocalGridType& getSubGrid(IndexType i) const;

public:
  /** Get the subgrid belonging to a given coarse cell.
  *
  * @param[in] coarseCell The coarse cell.
  * @return Returns the subgrid belonging to the given coarse cell.
  */
  const LocalGridType& getSubGrid(const CoarseEntityType& entity) const;
  LocalGridType& getSubGrid(const CoarseEntityType& entity);

  std::size_t size() const;

  //! returns true iff all corners of local_entity are inside coarse_entity
  static bool covers_strict(const CoarseEntityType& coarse_entity, const MsFEMTraits::LocalEntityType& local_entity);
  template <class GridImp, template <int, int, class> class GeometryImp>
  static bool covers_strict(
      const CoarseEntityType& coarse_entity,
      const Dune::Geometry<CommonTraits::world_dim, CommonTraits::world_dim, GridImp, GeometryImp>& local_geometry);
  //! returns true if local_entity's center is inside coarse_entity
  bool covers(const CoarseEntityType& coarse_entity, const MsFEMTraits::LocalEntityType& local_entity);

private:
  typedef std::map<IndexType, std::shared_ptr<LocalGridType>> SubGridStorageType;

  const CommonTraits::SpaceType& coarseSpace_;
  SubGridStorageType subGridList_;
  const LeafIndexSet& coarseGridLeafIndexSet_;
};

template <class GridImp, template <int, int, class> class GeometryImp>
bool LocalGridList::covers_strict(
    const CoarseEntityType& coarse_entity,
    const Dune::Geometry<CommonTraits::world_dim, CommonTraits::world_dim, GridImp, GeometryImp>& local_geometry) {
  const auto& reference_element = Stuff::Grid::reference_element(coarse_entity);
  const auto& coarse_geometry = coarse_entity.geometry();
  for (const auto i : DSC::valueRange(local_geometry.corners())) {
    if (!reference_element.checkInside(coarse_geometry.local(local_geometry.corner(i))))
      return false;
  }
  return true;
}


} // namespace Multiscale {
} // namespace Dune {

#endif // #ifndef SUBGRIDLIST_HH
