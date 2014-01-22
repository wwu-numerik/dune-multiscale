#include <config.h>
#include <assert.h>
#include <boost/assert.hpp>
#include <boost/multi_array/multi_array_ref.hpp>
#include <dune/common/exceptions.hh>
#include <dune/fem/misc/mpimanager.hh>
#include <dune/stuff/common/logging.hh>
#include <dune/stuff/common/profiler.hh>
#include <dune/stuff/common/ranges.hh>
#include <dune/stuff/common/float_cmp.hh>
#include <dune/stuff/grid/information.hh>
#include <algorithm>
#include <iterator>
#include <ostream>
#include <utility>
#include <memory>
#include <dune/multiscale/tools/misc.hh>
#include <dune/grid/common/gridfactory.hh>
#include <dune/stuff/grid/structuredgridfactory.hh>

#include <Eigen/Core>

#include "subgrid-list.hh"

namespace Dune {
namespace Multiscale {
namespace MsFEM {

LocalGridList::LocalGridList(MsFEMTraits::MacroMicroGridSpecifierType& specifier, bool silent /*= true*/)
  : coarseSpace_(specifier.coarseSpace())
  , specifier_(specifier)
  , silent_(silent)
  , coarseGridLeafIndexSet_(coarseSpace_.gridPart().grid().leafIndexSet())
  , fineToCoarseMap_(Fem::MPIManager::size()) {
  DSC::Profiler::ScopedTiming st("msfem.subgrid_list");
  createSubGrids();
}

LocalGridList::~LocalGridList() {}

/** Get the subgrid belonging to a given coarse cell index.
*
* @param[in] coarseCellIndex The index of a coarse cell.
* @return Returns the subgrid belonging to the coarse cell with the given index.
*/
MsFEMTraits::LocalGridType& LocalGridList::getSubGrid(std::size_t coarseCellIndex) {
  auto found = subGridList_.find(coarseCellIndex);
  BOOST_ASSERT_MSG(found != subGridList_.end(), "There is no subgrid for the index you provided!");
  assert(found->second);
  return *(found->second);
} // getSubGrid

/** Get the subgrid belonging to a given coarse cell index.
*
* @param[in] coarseCellIndex The index of a coarse cell.
* @return Returns the subgrid belonging to the coarse cell with the given index.
*/
const MsFEMTraits::LocalGridType& LocalGridList::getSubGrid(std::size_t coarseCellIndex) const {
  auto found = subGridList_.find(coarseCellIndex);
  BOOST_ASSERT_MSG(found != subGridList_.end(), "There is no subgrid for the index you provided!");
  assert(found->second);
  return *(found->second);
} // getSubGrid

/** Get the subgrid belonging to a given coarse cell.
*
* @param[in] coarseCell The coarse cell.
* @return Returns the subgrid belonging to the given coarse cell.
*/
const MsFEMTraits::LocalGridType& LocalGridList::getSubGrid(const CoarseEntityType& entity) const {
  const int index = coarseGridLeafIndexSet_.index(entity);
  return getSubGrid(index);
} // getSubGrid

/** Get the subgrid belonging to a given coarse cell.
*
* @param[in] coarseCell The coarse cell.
* @return Returns the subgrid belonging to the given coarse cell.
*/
MsFEMTraits::LocalGridType& LocalGridList::getSubGrid(const CoarseEntityType& entity) {
  const int index = coarseGridLeafIndexSet_.index(entity);
  return getSubGrid(index);
} // getSubGrid

// given the id of a subgrid, return the entity seed for the 'base coarse entity'
// (i.e. the coarse entity that the subgrid was constructed from by enrichment )
const LocalGridList::CoarseEntitySeedType &LocalGridList::get_coarse_entity_seed(std::size_t i) const {
  // the following returns the mapped element for index i if present,
  // if not, an out-of-range exception is thrown
  assert(false);//need to eliminate narrowing conversion
  return subgrid_id_to_base_coarse_entity_.at(i);
}


// only required for oversampling strategies with constraints (e.g strategy 2 or 3):
// for each given subgrid index return the vecor of ALL coarse nodes (global coordinates) that are in the subgrid,
// this also includes the coarse nodes on the boundary of U(T), even if this is a global Dirichlet node!
const LocalGridList::CoarseNodeVectorType& LocalGridList::getCoarseNodeVector(std::size_t i) const {
  if (DSC_CONFIG_GET("msfem.oversampling_strategy", 1) == 1)
    DUNE_THROW(Dune::InvalidStateException, "Method 'getCoarseNodeVector' of class 'LocalGridList' should not be used in\
                combination with oversampling strategy 1. Check your implementation!");

  if (i >= specifier_.getNumOfCoarseEntities()) {
    DUNE_THROW(Dune::RangeError, "Error. Subgrid-Index too large.");
  }
  return coarse_node_store_[i];
} // getSubGrid

// only required for oversampling strategy 3:
// this method only differs from the method 'getCoarseNodeVector' if the oversampling patch
// cannot be excluively described by a union of coarse grid elements.
// According to the definition of the LOD 'not full coarse layers' require that the averaging
// property of the weighted Clement operator is also applied to those coarse nodes, where
// the corresponding basis function has a nonempty intersection with the patch
const LocalGridList::CoarseNodeVectorType& LocalGridList::getExtendedCoarseNodeVector(std::size_t i) const {
  if ((DSC_CONFIG_GET("msfem.oversampling_strategy", 1) == 1) || (DSC_CONFIG_GET("msfem.oversampling_strategy", 1) == 2))
    DUNE_THROW(Dune::InvalidStateException,
               "Method 'getExtendendCoarseNodeVector' of class 'LocalGridList' should not be used in\
                combination with oversampling strategy 1 or 2. Check your implementation!");

  if (i >= specifier_.getNumOfCoarseEntities()) {
    DUNE_THROW(Dune::RangeError, "Error. Subgrid-Index too large.");
  }
  return extended_coarse_node_store_[i];
} // getSubGrid

// get number of sub grids
std::size_t LocalGridList::size() const { return specifier_.getNumOfCoarseEntities(); }

MsFEM::MsFEMTraits::LocalGridPartType LocalGridList::gridPart(std::size_t i) { return LocalGridPartType(getSubGrid(i)); }

bool LocalGridList::covers(const CoarseEntityType &coarse_entity, const LocalEntityType &local_entity) {
  const auto& center = local_entity.geometry().center();
  const auto& coarse_geo = coarse_entity.geometry();
  const auto center_local = coarse_geo.local(center);
  const auto& reference_element = Stuff::Grid::reference_element(coarse_entity);
  return reference_element.checkInside(center_local);
}

void LocalGridList::createSubGrids() {
//  DSC_PROFILER ("msfem.subgrid_list.create");
  DSC_LOG_INFO << "Starting creation of subgrids." << std::endl << std::endl;

  // the number of coarse grid entities (of codim 0).
  const auto number_of_coarse_grid_entities = specifier_.getNumOfCoarseEntities();
  const auto oversampling_strategy = DSC_CONFIG_GET("msfem.oversampling_strategy", 1);

  if ((oversampling_strategy == 2) || (oversampling_strategy == 3)) {
    coarse_node_store_ = CoarseGridNodeStorageType(number_of_coarse_grid_entities, CoarseNodeVectorType());
    extended_coarse_node_store_ = CoarseGridNodeStorageType(number_of_coarse_grid_entities, CoarseNodeVectorType());
  }
  typedef StructuredGridFactory<LocalGridType> FactoryType;
  const auto oversampling_layer = DSC_CONFIG_GET("msfem.oversampling_strategy", 1);
  const auto coarse_dimensions = DSG::dimensions<CommonTraits::GridType>(coarseSpace_.gridPart().grid());


  for (const auto& coarse_entity : coarseSpace_) {
    // make sure we only create subgrids for interior coarse elements, not
    // for overlap or ghost elements
    assert(coarse_entity.partitionType() == Dune::InteriorEntity);
    const auto coarse_index = coarseGridLeafIndexSet_.index(coarse_entity);
    // make sure we did not create a subgrid for the current coarse entity so far
    assert(subGridList_.find(coarse_index) == subGridList_.end());
    subgrid_id_to_base_coarse_entity_.insert(std::make_pair(coarse_index, std::move(coarse_entity.seed())));

    const auto dimensions = DSG::dimensions<CommonTraits::GridType>(coarse_entity);
    const int dim_world = LocalGridType::dimensionworld;
    typedef FieldVector<typename LocalGridType::ctype, dim_world> CoordType;
    CoordType lowerLeft(0);
    CoordType upperRight(0);
    array<unsigned int,dim_world> elemens;
    for(const auto i : DSC::valueRange(dim_world))
    {
      elemens[i] = DSC_CONFIG_GET("msfem.micro_cells_per_macrocell_dim", 8);
      const auto min = dimensions.coord_limits[i].min();
      const auto max = dimensions.coord_limits[i].max();
      const auto coarse_min = coarse_dimensions.coord_limits[i].min();
      const auto coarse_max = coarse_dimensions.coord_limits[i].max();
      const auto delta = max - min;
      lowerLeft[i] = std::max(min - oversampling_layer * delta, coarse_min);
      upperRight[i] = std::min(max + oversampling_layer * delta, coarse_max);
    }

    boost::format sp("Subgrid %d from (%f,%f) to (%f,%f) created.\n");
    DSC_LOG_INFO << sp % coarse_index % lowerLeft[0] % lowerLeft[1] % upperRight[0] % upperRight[1];
    subGridList_[coarse_index] = FactoryType::createCubeGrid(lowerLeft, upperRight, elemens);

    if ((oversampling_strategy == 2) || (oversampling_strategy == 3)) {
      assert(coarse_index >= 0 && coarse_index < coarse_node_store_.size() &&
             "Index set is not suitable for the current implementation!");
      for (int c = 0; c < coarse_entity.geometry().corners(); ++c) {
        coarse_node_store_[coarse_index].emplace_back(coarse_entity.geometry().corner(c));
        extended_coarse_node_store_[coarse_index].emplace_back(coarse_entity.geometry().corner(c));
      }
    }
  }
}

} // namespace MsFEM {
} // namespace Multiscale {
} // namespace Dune {
