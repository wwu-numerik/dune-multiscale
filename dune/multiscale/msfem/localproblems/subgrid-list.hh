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
#include <dune/fem/operator/2order/lagrangematrixsetup.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/space/lagrange.hh>
#include <dune/multiscale/msfem/msfem_grid_specifier.hh>
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
  typedef typename LocalGridType::Traits::LeafIndexSet LocalGridLeafIndexSet;
  typedef typename LocalGridLeafIndexSet::IndexType EntityIndexType;

  typedef typename LocalGridDiscreteFunctionSpaceType::IteratorType LocalGridEntityIteratorType;
  typedef typename LocalGridEntityIteratorType::Entity LocalEntityType;
  typedef typename LocalEntityType::EntityPointer LocalEntityPointerType;
  typedef typename LocalEntityType::EntitySeed LocalGridEntitySeed;
  typedef typename LocalEntityType::Codim<LocalGridType::dimension>::EntityPointer HostNodePointer;
  typedef typename LocalGridDiscreteFunctionSpaceType::GridPartType::IntersectionIteratorType HostIntersectionIterator;

  typedef typename MsFEMTraits::CoarseEntityType CoarseEntityType;
  typedef typename CoarseEntityType::EntitySeed CoarseEntitySeedType;
  typedef typename LeafIndexSet::IndexType CoarseEntityIndexType;
  typedef typename LocalGridDiscreteFunctionSpaceType::FunctionSpaceType::DomainType DomainType;

public:
  typedef typename LocalGridType::Traits::GlobalIdSet::IdType IdType;
  typedef std::vector<DomainType> CoarseNodeVectorType;

private:
  typedef std::vector<CoarseNodeVectorType> CoarseGridNodeStorageType; 
  typedef std::vector<std::vector<LocalEntityPointerType>> EntityPointerCollectionType;

public:
  LocalGridList(MsFEMTraits::MacroMicroGridSpecifierType& specifier, bool silent = true);
  ~LocalGridList();

private:
  LocalGridType& getSubGrid(IndexType i);
  const LocalGridType& getSubGrid(IndexType i) const;

public:
  const LocalGridType& getSubGrid(const CoarseEntityType& entity) const;
  LocalGridType& getSubGrid(const CoarseEntityType& entity);

  std::size_t size() const;

  LocalGridPartType gridPart(IndexType i);

  // given the index of a (codim 0) host grid entity, return the indices of the subgrids that contain the entity
  const std::vector<IndexType>& getSubgridIDs_that_contain_entity(IndexType host_enitity_index) const;

  // only required for oversampling strategies with constraints (e.g strategy 2 or 3):
  const CoarseNodeVectorType& getCoarseNodeVector(IndexType i) const;

  // only required for oversampling strategy 3:
  // this method only differs from the method 'getCoarseNodeVector' if the oversampling patch
  // cannot be excluively described by a union of coarse grid elements.
  // According to the definition of the LOD 'not full coarse layers' require that the averaging
  // property of the weighted Clement operator is also applied to those coarse nodes, where
  // the corresponding basis function has a nonempty intersection with the patch
  const CoarseNodeVectorType& getExtendedCoarseNodeVector(IndexType i) const;

  //! returns true if the coarse_entity covers the local_entity
  bool covers(const CoarseEntityType& coarse_entity, const LocalEntityType& local_entity);

  /** Get the mapping from node number to codim 0 host entity.
  * @return Returns the map.
  * */
  const EntityPointerCollectionType& getNodeEntityMap();

  // given the id of a subgrid, return the entity seed for the 'base coarse entity'
  // (i.e. the coarse entity that the subgrid was constructed from by enrichment )
  const CoarseEntitySeedType& get_coarse_entity_seed(std::size_t i) const;

private:
  typedef std::map<IndexType, std::shared_ptr<LocalGridType>> SubGridStorageType;

  void createSubGrids();

  const CommonTraits::DiscreteFunctionSpaceType& coarseSpace_;
  MsFEMTraits::MacroMicroGridSpecifierType& specifier_;
  bool silent_;
  SubGridStorageType subGridList_;
  CoarseGridNodeStorageType coarse_node_store_;
  CoarseGridNodeStorageType extended_coarse_node_store_;
  const LeafIndexSet& coarseGridLeafIndexSet_;
  std::vector<std::map<IndexType, IndexType>> fineToCoarseMap_;
  std::map<IdType, IdType> fineToCoarseMapID_;

  // given the id of a subgrid, return the entity seed for the 'base coarse entity'
  // (i.e. the coarse entity that the subgrid was constructed from by enrichment )
  std::map<CoarseEntityIndexType, CoarseEntitySeedType> subgrid_id_to_base_coarse_entity_;
};

} // namespace MsFEM {
} // namespace Multiscale {
} // namespace Dune {

#endif // #ifndef SUBGRIDLIST_HH
