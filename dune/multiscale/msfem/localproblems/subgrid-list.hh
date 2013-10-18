// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef SUBGRIDLIST_HH
#define SUBGRIDLIST_HH


#include <boost/noncopyable.hpp>
#include <boost/multi_array.hpp>

#include <dune/common/fmatrix.hh>
#include <dune/common/shared_ptr.hh>
#include <dune/common/exceptions.hh>
#include <dune/stuff/grid/entity.hh>

#include <dune/multiscale/tools/subgrid_io.hh>
#include <dune/subgrid/subgrid.hh>

#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/space/lagrange.hh>
//#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/operator/2order/lagrangematrixsetup.hh>

#include <dune/multiscale/msfem/msfem_traits.hh>
#include <dune/multiscale/msfem/msfem_grid_specifier.hh>

namespace Dune {
namespace Multiscale {
namespace MsFEM {

//! container for cell problem subgrids
class SubGridList : public boost::noncopyable {
  typedef typename CommonTraits::DiscreteFunctionType HostDiscreteFunctionImp;

public:
  //! ---------------- typedefs for the HostDiscreteFunctionSpace -----------------------

  typedef MsFEMTraits::MacroMicroGridSpecifierType MacroMicroGridSpecifierType;
  typedef HostDiscreteFunctionImp HostDiscreteFunctionType;
  //! type of discrete function space
  typedef typename HostDiscreteFunctionType::DiscreteFunctionSpaceType HostDiscreteFunctionSpaceType;
  //! type of grid partition
  typedef typename HostDiscreteFunctionSpaceType::GridPartType HostGridPartType;

  //! type of grid
private:
  typedef typename HostDiscreteFunctionSpaceType::GridType HostGridType;
  typedef typename HostGridType::Traits::LeafIndexSet HostGridLeafIndexSet;
  typedef typename HostGridLeafIndexSet::IndexType EntityIndexType;
  typedef typename HostDiscreteFunctionSpaceType::IteratorType HostGridEntityIteratorType;
  typedef typename HostGridEntityIteratorType::Entity HostEntityType;
  typedef typename HostEntityType::EntityPointer HostEntityPointerType;
  typedef typename HostEntityType::EntitySeed CoarseGridEntitySeed;
  typedef typename HostEntityType::Codim<HostGridType::dimension>::EntityPointer HostNodePointer;
  typedef typename HostGridPartType::IntersectionIteratorType HostIntersectionIterator;

  typedef typename MsFEMTraits::CoarseEntityType CoarseEntityType;

  //! type of (non-discrete )function space
  typedef typename HostDiscreteFunctionSpaceType::FunctionSpaceType FunctionSpaceType;

  //! type of domain
  typedef typename FunctionSpaceType::DomainType DomainType;

public:
  typedef typename HostGridType::Traits::GlobalIdSet::IdType IdType;
  typedef std::vector<DomainType> CoarseNodeVectorType;

private:
  typedef std::vector<CoarseNodeVectorType> CoarseGridNodeStorageType;
  typedef boost::multi_array<bool, 3> EnrichmentMatrixType;
  //! @todo this should eventually be changed to the type of the coarse space
  typedef typename HostGridPartType::Codim<0>::EntityType::Geometry::LocalCoordinate LocalCoordinateType;
  typedef ReferenceElements<typename LocalCoordinateType::value_type, LocalCoordinateType::dimension>
  CoarseRefElementType;
  typedef std::vector<std::vector<HostEntityPointerType>> EntityPointerCollectionType;

public:
  //! ---------------- typedefs for the SubgridDiscreteFunctionSpace -----------------------
  // ( typedefs for the local grid and the corresponding local ('sub') )discrete space )
  //! type of grid
  typedef MsFEMTraits::SubGridType SubGridType;
  //! type of grid part
  typedef MsFEMTraits::SubGridPartType SubGridPartType;
  //! type of subgrid discrete function space
  typedef MsFEMTraits::SubGridDiscreteFunctionSpaceType SubGridDiscreteFunctionSpaceType;
  //! type of subgrid discrete function
  typedef MsFEMTraits::SubGridDiscreteFunctionType SubGridDiscreteFunctionType;

  SubGridList(MacroMicroGridSpecifierType& specifier, bool silent = true);
  ~SubGridList();

private:
  SubGridType& getSubGrid(std::size_t i);
  const SubGridType& getSubGrid(std::size_t i) const;

public:
  const SubGridType& getSubGrid(const CoarseEntityType& entity) const;
  SubGridType& getSubGrid(const CoarseEntityType& entity);

  std::size_t DUNE_DEPRECATED_MSG("use size() instead") getNumberOfSubGrids() const;
  std::size_t size() const;

  SubGridPartType gridPart(std::size_t i);

  // given the index of a (codim 0) host grid entity, return the indices of the subgrids that contain the entity
  const std::vector<std::size_t>& getSubgridIDs_that_contain_entity(std::size_t host_enitity_index) const;

  // only required for oversampling strategies with constraints (e.g strategy 2 or 3):
  const CoarseNodeVectorType& getCoarseNodeVector(std::size_t i) const;

  // only required for oversampling strategy 3:
  // this method only differs from the method 'getCoarseNodeVector' if the oversampling patch
  // cannot be excluively described by a union of coarse grid elements.
  // According to the definition of the LOD 'not full coarse layers' require that the averaging
  // property of the weighted Clement operator is also applied to those coarse nodes, where
  // the corresponding basis function has a nonempty intersection with the patch
  const CoarseNodeVectorType& getExtendedCoarseNodeVector(std::size_t i) const;

  /** Get the index of the coarse cell enclosing the barycentre of a given fine cell.
*
* Given a fine cell, this method computes its barycentre. Using a grid run on the coarse
* grid, it checks which (if any) coarse cell contains the barycentre.
*
* @param[in] hostEntity The host entity.
* @param[in,out] lastIterator The macro cell that was found in the last run. This should be set to
*                             coarseGrid.begin<0>() in the first run. This iterator will then be
*                             updated and set to the macro element used in this run.
* @param[in] coarseGridLeafIndexSet The index set of the coarse grid.
*
*/
  std::size_t getEnclosingMacroCellIndex(const HostEntityPointerType& hostEntityPointer);

  IdType getEnclosingMacroCellId(const HostEntityPointerType& hostEntityPointer);

  /** Get the mapping from node number to codim 0 host entity.
  * @return Returns the map.
  * */
  const EntityPointerCollectionType& getNodeEntityMap();

  // given the id of a subgrid, return the entity seed for the 'base coarse entity'
  // (i.e. the coarse entity that the subgrid was constructed from by enrichment )
  const CoarseGridEntitySeed& get_coarse_entity_seed(std::size_t i) const;

private:
  typedef std::map<std::size_t, std::shared_ptr<SubGridType>> SubGridStorageType;
  /**
   * \note called in SubGridList constructor only
   */
  void enrichment(const HostEntityPointerType& hit,
                  //          const HostEntityPointerType& level_father_it,
                  const std::size_t& father_index, // father_index = index/number of current subgrid
                  std::shared_ptr<SubGridType> subGrid, std::size_t& layer);
  bool entityPatchInSubgrid(const HostEntityPointerType& hit, const HostGridPartType& hostGridPart,
                            shared_ptr<const SubGridType> subGrid,
                            const EntityPointerCollectionType& entities_sharing_same_node) const;

  void identifySubGrids();
  void createSubGrids();
  void finalizeSubGrids();

  const HostDiscreteFunctionSpaceType& hostSpace_;
  const HostDiscreteFunctionSpaceType& coarseSpace_;
  MacroMicroGridSpecifierType& specifier_;
  bool silent_;
  SubGridStorageType subGridList_;
  CoarseGridNodeStorageType coarse_node_store_;
  CoarseGridNodeStorageType extended_coarse_node_store_;
  const HostGridLeafIndexSet& coarseGridLeafIndexSet_;
  const HostGridLeafIndexSet& hostGridLeafIndexSet_;
  const HostGridPartType& hostGridPart_;
  EntityPointerCollectionType entities_sharing_same_node_;
  EnrichmentMatrixType enriched_;
  std::vector<std::map<std::size_t, std::size_t>> fineToCoarseMap_;
  std::map<IdType, IdType> fineToCoarseMapID_;
  // given the id of a fine grid element, the vector returns the ids of all subgrids that share that element
  std::vector<std::vector<std::size_t>> fine_id_to_subgrid_ids_;

  // given the id of a subgrid, return the entity seed for the 'base coarse entity'
  // (i.e. the coarse entity that the subgrid was constructed from by enrichment )
  std::map<EntityIndexType, CoarseGridEntitySeed> subgrid_id_to_base_coarse_entity_;
};

} // namespace MsFEM {
} // namespace Multiscale {
} // namespace Dune {

#endif // #ifndef SUBGRIDLIST_HH
