#ifndef DUNE_FEM_PERIODICGRIDPART_ITERATOR_HH
#define DUNE_FEM_PERIODICGRIDPART_ITERATOR_HH


template< class Grid >
class PeriodicLeafIntersectionIterator
{
  typedef PeriodicLeafIntersectionIterator< Grid > ThisType;

  friend class PeriodicLeafGridPart< Grid >;

public:
  typedef Grid GridType;

  typedef ThisType Intersection;

protected:
  typedef typename GridType::template Codim< 0 >::Entity Codim0EntityType;

  typedef typename Codim0EntityType::LeafIntersectionIterator
  WrappedIteratorType;

public:
  typedef typename WrappedIteratorType::ctype ctype;

  enum
  {
    dimension = WrappedIteratorType::dimension,
    dimensionworld = WrappedIteratorType::dimensionworld
  };

  typedef FieldVector< ctype, dimensionworld > DomainType;

  typedef FieldVector< ctype, dimensionworld - 1 > LocalDomainType;

  typedef typename WrappedIteratorType::Entity Entity;

  typedef typename WrappedIteratorType::EntityPointer EntityPointer;

  typedef typename WrappedIteratorType::Geometry Geometry;

  typedef typename WrappedIteratorType::LocalGeometry LocalGeometry;

  typedef typename WrappedIteratorType::ImplementationType
  ImplementationType;

protected:
  WrappedIteratorType wrappedIterator_;

protected:
  explicit PeriodicLeafIntersectionIterator(WrappedIteratorType wrappedIterator)
    : wrappedIterator_(wrappedIterator)
  {}

public:
  PeriodicLeafIntersectionIterator(const ThisType& other)
    : wrappedIterator_(other.wrappedIterator_)
  {}

  ThisType& operator=(const ThisType& other) {
    wrappedIterator_ = other.wrappedIterator_;
    return *this;
  }

  ThisType& operator++() {
    ++wrappedIterator_;
    return *this;
  }

  const Intersection& operator*() const {
    return *this;
  }

  const Intersection* operator->() const {
    return this;
  }

  bool operator==(const ThisType& other) const {
    return wrappedIterator_ == other.wrappedIterator_;
  }

  bool operator!=(const ThisType& other) const {
    return wrappedIterator_ != other.wrappedIterator_;
  }

  bool boundary() const {
    // IntersectionIterator specifies that this should return true at the boundary
    return false;
  }

  int boundaryId() const {
    return 0;
  }

  int neighbor() const {
    // IntersectionIterator specifies that we should return true!
    return wrappedIterator_->neighbor();
  }

  EntityPointer inside() const {
    return wrappedIterator_->inside();
  }

  EntityPointer outside() const {
    if ( wrappedIterator_->neighbor() )
      return wrappedIterator_->outside();
    else
      DUNE_THROW(NotImplemented, "PeriodicLeafIntersectionIteratorWrapper: "
                                 "outside on boundary not implemented yet.");
  } // outside

  const LocalGeometry& intersectionSelfLocal() const {
    return wrappedIterator_->intersectionSelfLocal();
  }

  const LocalGeometry& intersectionNeighborLocal() const {
    if ( wrappedIterator_->neighbor() )
      return wrappedIterator_->intersectionNeighborLocal();
    else
      DUNE_THROW(NotImplemented, "PeriodicLeafIntersectionIteratorWrapper: "
                                 "outside on boundary not implemented yet.");
  } // intersectionNeighborLocal

  const Geometry& intersectionGlobal() const {
    return wrappedIterator_->intersectionGlobal();
  }

  const LocalGeometry& geometryInInside() const {
    return wrappedIterator_->geometryInInside();
  }

  const LocalGeometry& geometryInOutside() const {
    if ( wrappedIterator_->neighbor() )
      return wrappedIterator_->geometryInOutside();
    else
      DUNE_THROW(NotImplemented, "PeriodicLeafIntersectionIteratorWrapper: "
                                 "outside on boundary not implemented yet.");
  } // geometryInOutside

  const Geometry& geometry() const {
    return wrappedIterator_->geometry();
  }

  int numberInSelf() const {
    return wrappedIterator_->numberInSelf();
  }

  int numberInNeighbor() const {
    return wrappedIterator_->numberInNeighbor();
  }

  int indexInInside() const {
    return wrappedIterator_->indexInInside();
  }

  int indexInOutside() const {
    return wrappedIterator_->indexInOutside();
  }

  DomainType outerNormal(const LocalDomainType& x) const {
    return wrappedIterator_->outerNormal(x);
  }

  DomainType integrationOuterNormal(const LocalDomainType& x) const {
    return wrappedIterator_->integrationOuterNormal(x);
  }

  DomainType unitOuterNormal(const LocalDomainType& x) const {
    return wrappedIterator_->unitOuterNormal(x);
  }

protected:
  const ImplementationType& getRealImp() const {
    return wrappedIterator_.getRealImp();
  }

  ImplementationType& getRealImp() {
    return wrappedIterator_.getRealImp();
  }
};


#endif // DUNE_FEM_PERIODICGRIDPART_ITERATOR_HH
