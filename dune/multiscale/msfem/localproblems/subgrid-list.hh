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
#include <dune/fem/space/lagrangespace.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/operator/2order/lagrangematrixsetup.hh>

#include <dune/multiscale/msfem/msfem_traits.hh>
#include <dune/multiscale/msfem/msfem_grid_specifier.hh>

namespace Dune {
namespace Multiscale {
namespace MsFEM {

//! container for cell problem subgrids
class SubGridList : public boost::noncopyable
{
  typedef typename CommonTraits::DiscreteFunctionType HostDiscreteFunctionImp;
  typedef MsFEMTraits::SubGridType SubGridImp;
  typedef MsFEMTraits::MacroMicroGridSpecifierType MacroMicroGridSpecifierImp;

public:
  //! ---------------- typedefs for the HostDiscreteFunctionSpace -----------------------

  typedef MacroMicroGridSpecifierImp MacroMicroGridSpecifierType;
  typedef HostDiscreteFunctionImp HostDiscreteFunctionType;
  //! type of discrete function space
  typedef typename HostDiscreteFunctionType::DiscreteFunctionSpaceType HostDiscreteFunctionSpaceType;
  //! type of grid partition
  typedef typename HostDiscreteFunctionSpaceType::GridPartType HostGridPartType;

  //! type of grid
private:
  typedef typename HostDiscreteFunctionSpaceType::GridType HostGridType;
  typedef typename HostGridType::Traits::LeafIndexSet HostGridLeafIndexSet;
  typedef typename HostDiscreteFunctionSpaceType::IteratorType HostGridEntityIteratorType;
  typedef typename HostGridEntityIteratorType::Entity HostEntityType;
  typedef typename HostEntityType::EntityPointer HostEntityPointerType;
  typedef typename HostEntityType::Codim< 2 >::EntityPointer HostNodePointer;
  typedef typename HostGridPartType::IntersectionIteratorType HostIntersectionIterator;

  //! type of (non-discrete )function space
  typedef typename HostDiscreteFunctionSpaceType::FunctionSpaceType FunctionSpaceType;

  //! type of domain
  typedef typename FunctionSpaceType::DomainType DomainType;

public:
  typedef std::vector< DomainType > CoarseNodeVectorType;

private:
  typedef std::vector< CoarseNodeVectorType > CoarseGridNodeStorageType;
  typedef boost::multi_array<bool, 3> EnrichmentMatrixType;
  //! @todo this should eventually be changed to the type of the coarse space
  typedef typename HostGridPartType::Codim<0>::EntityType::Geometry::LocalCoordinate LocalCoordinateType;
  typedef GenericReferenceElements< typename  LocalCoordinateType::value_type,  LocalCoordinateType::dimension >
          CoarseRefElementType;
  typedef std::vector<std::vector<HostEntityPointerType>> EntityPointerCollectionType;

public:
  //! ---------------- typedefs for the SubgridDiscreteFunctionSpace -----------------------
  // ( typedefs for the local grid and the corresponding local ('sub') )discrete space )

  //! type of grid
  typedef SubGridImp SubGridType;

  //! type of grid part
  typedef LeafGridPart< SubGridType > SubGridPartType;
  
    //! type of subgrid discrete function space
  typedef LagrangeDiscreteFunctionSpace< FunctionSpaceType, SubGridPartType, 1/*=POLORDER*/ > SubGridDiscreteFunctionSpace;

  //! type of subgrid discrete function
  typedef AdaptiveDiscreteFunction< SubGridDiscreteFunctionSpace > SubGridDiscreteFunction;

  SubGridList(MacroMicroGridSpecifierType& specifier, bool silent = true);
  ~SubGridList();
  SubGridType& getSubGrid(int i);
  const SubGridType& getSubGrid(int i) const;

  // only required for oversampling strategies with constraints (e.g strategy 2 or 3):
  const CoarseNodeVectorType& getCoarseNodeVector(int i) const;

private:
  typedef std::vector< std::shared_ptr<SubGridType> > SubGridStorageType;
  /**
   * \note called in SubGridList constructor only
   */
  void enrichment(const HostEntityPointerType& hit,
          const HostEntityPointerType& level_father_it,
          const int& father_index, // father_index = index/number of current subgrid
          shared_ptr<SubGridType> subGrid,
          int& layer);

  bool entityPatchInSubgrid(const HostEntityPointerType& hit,
          const HostGridPartType& hostGridPart,
          shared_ptr<const SubGridType> subGrid,
          const EntityPointerCollectionType& entities_sharing_same_node) const;

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
  int getEnclosingMacroCellIndex(const HostEntityPointerType& hostEntityPointer);

  void identifySubGrids();
  void createSubGrids();
  void finalizeSubGrids();

  const HostDiscreteFunctionSpaceType& hostSpace_;
  const HostDiscreteFunctionSpaceType& coarseSpace_;
  MacroMicroGridSpecifierType& specifier_;
  bool silent_;
  SubGridStorageType subGridList_;
  CoarseGridNodeStorageType coarse_node_store_;
  const HostGridLeafIndexSet& coarseGridLeafIndexSet_;
  const HostGridLeafIndexSet& hostGridLeafIndexSet_;
  const HostGridPartType& hostGridPart_;
  EntityPointerCollectionType entities_sharing_same_node_;
  EnrichmentMatrixType enriched_;
  std::map<int, int> fineToCoarseMap_;
};

} //namespace MsFEM {
} //namespace Multiscale {
} //namespace Dune {

#endif // #ifndef SUBGRIDLIST_HH

