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
  typedef typename HostEntityType::template Codim< 2 >::EntityPointer HostNodePointer;
  typedef typename HostGridPartType::IntersectionIteratorType HostIntersectionIterator;

  //! type of (non-discrete )function space
  typedef typename HostDiscreteFunctionSpaceType::FunctionSpaceType FunctionSpaceType;

  //! type of domain
  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef std::vector< DomainType > CoarseNodeVectorType;
  typedef std::vector< CoarseNodeVectorType > CoarseGridNodeStorageType;
  typedef boost::multi_array<bool, 3> EnrichmentMatrixType;
  //! @todo this should eventually be changed to the type of the coarse space
  typedef typename HostGridPartType::template Codim<0>::EntityType::Geometry::LocalCoordinate LocalCoordinateType;
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
  SubGridType& getSubGrid(int i);
  const SubGridType& getSubGrid(int i) const;

  // only required for oversampling strategies with constraints (e.g strategy 2 or 3):
  const CoarseNodeVectorType& getCoarseNodeVector(int i) const;

private:
  /**
   * \note called in SubGridList constructor only
   */
  void enrichment(const HostEntityPointerType& hit,
          const HostEntityPointerType& level_father_it,
          const int& father_index, // father_index = index/number of current subgrid
          const HostGridPartType& hostGridPart,
          shared_ptr<SubGridType> subGrid,
          EntityPointerCollectionType& entities_sharing_same_node,
          int& layer,
          EnrichmentMatrixType& enriched);

  bool entity_patch_in_subgrid(const HostEntityPointerType& hit,
          const HostGridPartType& hostGridPart,
          shared_ptr<const SubGridType> subGrid,
          const EntityPointerCollectionType& entities_sharing_same_node) const;

  const HostDiscreteFunctionSpaceType& hostSpace_;
  MacroMicroGridSpecifierType& specifier_;

  bool silent_;

  typedef std::vector< std::shared_ptr<SubGridType> > SubGridStorageType;
  SubGridStorageType subGridList_;

  CoarseGridNodeStorageType coarse_node_store_;  

public:
  ~SubGridList()
  {}
};

} //namespace MsFEM {
} //namespace Multiscale {
} //namespace Dune {

#endif // #ifndef SUBGRIDLIST_HH

