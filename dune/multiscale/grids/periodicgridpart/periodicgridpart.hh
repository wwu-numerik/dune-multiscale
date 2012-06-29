#ifndef DUNE_FEM_PERIODICGRIDPART_HH
#define DUNE_FEM_PERIODICGRIDPART_HH

#include <dune/grid/common/grid.hh>

#include <dune/fem/gridpart/common/gridpart.hh>
#include "periodicindexset.hh"

namespace Dune {
template< class Grid >
class PeriodicLeafGridPart;

template< class Grid >
class PeriodicLeafGridPartTraits
{
  typedef PeriodicLeafGridPartTraits< Grid > ThisType;

public:
  typedef Grid GridType;

  typedef PeriodicLeafGridPart< GridType > GridPartType;

  typedef PeriodicLeafIndexSet< GridType > IndexSetType;

  static const PartitionIteratorType indexSetPartitionType = All_Partition;

  // typedef PeriodicLeafIntersectionIterator< GridType > IntersectionIteratorType;
  typedef typename GridType::LeafIntersectionIterator IntersectionIteratorType;

  template< int codim >
  struct Codim
  {
    // type of the entity (should be wrapped, too)
    typedef typename GridType::template Codim< codim >::Entity EntityType;

    typedef typename GridType::template Codim< codim >::Geometry      GeometryType;
    typedef typename GridType::template Codim< codim >::LocalGeometry LocalGeometryType;

    typedef typename GridType::template Codim< codim >::EntityPointer EntityPointerType;

    typedef typename GridType::template Codim< codim >::EntitySeed EntitySeedType;

    template< PartitionIteratorType pitype >
    struct Partition
    {
      typedef typename GridType
        ::template Codim< codim >::template Partition< pitype >::LeafIterator
      IteratorType;
    };
  };

  // ! is the grid partition conforming?
  static const bool conforming = Capabilities::isLeafwiseConforming< GridType >::v;
};

/*! \addtogroup PeriodicGridPart
   *  \class PeriodicLeafGridPart
   *  \brief A grid partition identifying opposite faces of the unit cube
   *
   *  Using any grid of the unit cube, this grid partition makes the grid
   *  periodic by identifying the indices of subentities of opposite faces.
   *
   *  \note Since the underlying grid does not know about the periodic boundary,
   *        refinement may break conformity (global refinement should work,
   *        though).
   *
   *  \note This grid partition says that there is no boundary. In DUNE, however,
   *        periodic boundaries shall be implemented as boundaries with ghost
   *        entities (however the FEM codes usually only check if an intersection
   *        is a boundary.
   *
   *  \todo The entity needs also to be wrapped, so that hasBoundaryIntersections
   *        always returns false
   *
   *  \todo Return correct neighbors for entities with boundary intersections.
   *
   *  \newimplementation Allows to construct globally refined grids for the
   *                     unitcube with periodic boundaries.
   */
template< class Grid >
class PeriodicLeafGridPart
  : public GridPartDefault< PeriodicLeafGridPartTraits< Grid > >
{
  typedef PeriodicLeafGridPart< Grid >                          ThisType;
  typedef GridPartDefault< PeriodicLeafGridPartTraits< Grid > > BaseType;

public:
  // ! type of traits
  typedef typename BaseType::Traits Traits;

  // ! type of the underlying grid
  typedef typename Traits::GridType GridType;

    typedef typename GridType :: template Partition < All_Partition > :: LeafGridView LeafGridView;

  // ! type of the index set
  typedef typename Traits::IndexSetType IndexSetType;

  // ! type of intersection iterators
  typedef typename Traits::IntersectionIteratorType IntersectionIteratorType;

  template< int codim >
  struct Codim
    : public BaseType::template Codim< codim >
  {
    typedef typename Traits::template Codim< codim >::EntityType EntityType;
  };

  // ! Constructor wrapping a grid in the grid partition
  PeriodicLeafGridPart(GridType& grid)
    : BaseType(grid)
      , indexSet_(grid)
  {}

  // ! Returns reference to index set of the underlying grid
  const IndexSetType& indexSet() const {
    return indexSet_;
  }

  // ! Begin iterator on the leaf level
  template< int codim >
  typename Codim< codim >::IteratorType begin() const {
    return begin< codim, InteriorBorder_Partition >();
  }

  // ! Begin iterator on the leaf level
  template< int codim, PartitionIteratorType pitype >
  typename Codim< codim >::template Partition< pitype >::IteratorType begin() const {
    return (*this).grid().template leafbegin< codim, pitype >();
  }

  // ! Begin iterator on the leaf level
  template< int codim >
  typename Codim< codim >::IteratorType end() const {
    return end< codim, InteriorBorder_Partition >();
  }

  // ! End iterator on the leaf level
  template< int codim, PartitionIteratorType pitype >
  typename Codim< codim >::template Partition< pitype >::IteratorType end() const {
    return (*this).grid().template leafend< codim, pitype >();
  }

  /** \brief begin intersection iterator for an entity
     *
     *  \note The intersection iterators always return boundary = false
     *        and neighbor = true.
     *
     *  \param[in]  entity  entity the intersection iterator is requested for
     *
     *  \returns a begin intersection iterator
     */
  IntersectionIteratorType ibegin(const typename Codim< 0 >::EntityType& entity) const {
    return IntersectionIteratorType( entity.ileafbegin() );
  }

  /** \brief end intersection iterator for an entity
     *
     *  \note The intersection iterators always return boundary = false
     *        and neighbor = true.
     *
     *  \param[in]  entity  entity the intersection iterator is requested for
     *
     *  \returns an end intersection iterator
     */
  IntersectionIteratorType iend(const typename Codim< 0 >::EntityType& entity) const {
    return IntersectionIteratorType( entity.ileafend() );
  }

  // ! Deliver maximum level of grid
  int level() const {
    return this->grid().maxLevel();
  }

  // ! Communication Method for this grid partition
  template< class DataHandleImp, class DataType >
  void communicate(CommDataHandleIF< DataHandleImp, DataType >& data,
                   InterfaceType iftype,
                   CommunicationDirection dir) const {
    this->grid().communicate(data, iftype, dir);
  }

protected:
  IndexSetType indexSet_;
  // the leaf grid view
//  LeafGridView leafView_ ;
};
}

#endif // ifndef DUNE_FEM_PERIODICGRIDPART_HH
