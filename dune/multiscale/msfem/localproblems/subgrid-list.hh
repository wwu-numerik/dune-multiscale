// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef SUBGRIDLIST_HH
#define SUBGRIDLIST_HH

#include <boost/noncopyable.hpp>
#include <dune/common/deprecated.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/shared_ptr.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
//#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/space/lagrange.hh>
#include <dune/multiscale/msfem/msfem_traits.hh>
#include <dune/stuff/grid/entity.hh>

#include <cstddef>
#include <map>
#include <memory>
#include <vector>

#include "dune/multiscale/common/traits.hh"

namespace Dune {
namespace Multiscale {
namespace MsFEM {

//! container for cell problem subgrids
class LocalGridList : public boost::noncopyable {
  typedef typename MsFEMTraits::LocalGridDiscreteFunctionType LocalGridDiscreteFunctionType;
  typedef typename LocalGridDiscreteFunctionType::DiscreteFunctionSpaceType LocalGridDiscreteFunctionSpaceType;
  typedef typename LocalGridDiscreteFunctionSpaceType::GridType LocalGridType;
  typedef typename LocalGridDiscreteFunctionSpaceType::GridPartType LocalGridPartType;
  typedef typename CommonTraits::GridType::Traits::LeafIndexSet LeafIndexSet;
  typedef typename LeafIndexSet::IndexType IndexType;

  typedef typename MsFEMTraits::CoarseEntityType CoarseEntityType;
  typedef typename CoarseEntityType::EntitySeed CoarseEntitySeedType;
  typedef typename LeafIndexSet::IndexType CoarseEntityIndexType;
  typedef typename LocalGridDiscreteFunctionSpaceType::FunctionSpaceType::DomainType DomainType;

public:
  typedef std::vector<DomainType> CoarseNodeVectorType;

  LocalGridList(const CommonTraits::DiscreteFunctionSpaceType& coarseSpace);

private:
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

  LocalGridPartType gridPart(IndexType i);
  //! returns true iff all corners of local_entity are inside coarse_entity
  bool covers_strict(const CoarseEntityType& coarse_entity, const MsFEMTraits::LocalEntityType& local_entity);
  template <class PointIterator>
  bool covers_strict(const CoarseEntityType& coarse_entity, const PointIterator first, const PointIterator last);
  //! returns true if local_entity's center is inside coarse_entity
  bool covers(const CoarseEntityType& coarse_entity, const MsFEMTraits::LocalEntityType& local_entity);

private:
  typedef std::map<IndexType, std::shared_ptr<LocalGridType>> SubGridStorageType;

  const CommonTraits::DiscreteFunctionSpaceType& coarseSpace_;
  SubGridStorageType subGridList_;
  const LeafIndexSet& coarseGridLeafIndexSet_;
};

} // namespace MsFEM {
} // namespace Multiscale {
} // namespace Dune {

#endif // #ifndef SUBGRIDLIST_HH
