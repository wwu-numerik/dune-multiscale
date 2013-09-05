#include "dirichletconstraints.hh"

namespace Dune {
namespace Multiscale {

template < class DomainSpace, class RangeSpace >
DirichletConstraints<DomainSpace, RangeSpace>::DirichletConstraints(const BoundaryType &boundary,
                                                                    const DomainSpaceType& domain_space )
: boundary_(boundary),
domain_space_( domain_space ),
slaveDofs_( getSlaveDofs( domain_space_ ) ),
dirichletBlocks_(),
// mark DoFs on the Dirichlet boundary
hasDirichletDofs_( false ),
sequence_( -1 )
{}

template < class DomainSpace, class RangeSpace >
template < class DiscreteFunctionType >
void DirichletConstraints<DomainSpace, RangeSpace>::operator()(const DiscreteFunctionType& u,
                                                               DiscreteFunctionType& w ) const
{
  updateDirichletDofs();

  // if Dirichlet Dofs have been found, treat them
  if( hasDirichletDofs_ )
  {
    typedef typename DiscreteFunctionType :: DofIteratorType DofIteratorType ;
    typedef typename DiscreteFunctionType :: ConstDofIteratorType ConstDofIteratorType ;

    ConstDofIteratorType uIt = u.dbegin();
    DofIteratorType wIt = w.dbegin();

    const unsigned int localBlockSize = DiscreteFunctionType :: DiscreteFunctionSpaceType ::
    localBlockSize ;
    // loop over all blocks
    const unsigned int blocks = u.space().blockMapper().size();
    for( unsigned int blockDof = 0; blockDof < blocks ; ++ blockDof )
    {
      if( dirichletBlocks_[ blockDof ] )
      {
        // copy dofs of the block
        for( unsigned int l = 0; l < localBlockSize ; ++ l, ++ wIt, ++ uIt )
        {
          assert( uIt != u.dend() );
          assert( wIt != w.dend() );
          (*wIt) = (*uIt);
        }
      }
      else
      {
        // increase dof iterators anyway
        for( unsigned int l = 0; l < localBlockSize ; ++ l, ++ wIt, ++ uIt )
        {}
      }
    }
  }
}


template < class DomainSpace, class RangeSpace >
template < class DiscreteFunctionType >
void DirichletConstraints<DomainSpace, RangeSpace>::setValue(const typename DiscreteFunctionType::RangeFieldType val,
                                                             DiscreteFunctionType& w ) const
{
  updateDirichletDofs();

  // if Dirichlet Dofs have been found, treat them
  if( hasDirichletDofs_ )
  {
    typedef typename DiscreteFunctionType :: DofIteratorType DofIteratorType ;
    typedef typename DiscreteFunctionType :: ConstDofIteratorType ConstDofIteratorType ;

    DofIteratorType wIt = w.dbegin();

    const unsigned int localBlockSize = DiscreteFunctionType :: DiscreteFunctionSpaceType ::
    localBlockSize ;
    // loop over all blocks
    const unsigned int blocks = w.space().blockMapper().size();
    for( unsigned int blockDof = 0; blockDof < blocks ; ++ blockDof )
    {
      if( dirichletBlocks_[ blockDof ] )
      {
        // copy dofs of the block
        for( unsigned int l = 0; l < localBlockSize ; ++ l, ++ wIt )
        {
          assert( wIt != w.dend() );
          (*wIt) = val;
        }
      }
      else
      {
        // increase dof iterators anyway
        for( unsigned int l = 0; l < localBlockSize ; ++ l, ++ wIt )
        {}
      }
    }
  }
}

template < class DomainSpace, class RangeSpace >
template < class GridFunctionType, class DiscreteFunctionType >
void DirichletConstraints<DomainSpace, RangeSpace>::operator ()(const GridFunctionType& u,
                                                                DiscreteFunctionType& w ) const
{
  apply( u, w );
}

template < class DomainSpace, class RangeSpace >
template <class LinearOperator>
void DirichletConstraints<DomainSpace, RangeSpace>::applyToOperator( LinearOperator& linearOperator ) const
{
  updateDirichletDofs();

  typedef typename DomainSpaceType :: IteratorType IteratorType;
  typedef typename IteratorType :: Entity EntityType;

  // if Dirichlet Dofs have been found, treat them
  if( hasDirichletDofs_ )
  {
    const IteratorType end = domain_space_.end();
    for( IteratorType it = domain_space_.begin(); it != end; ++it )
    {
      const EntityType &entity = *it;
      // adjust linear operator
      dirichletDofsCorrectOnEntity( linearOperator, entity );
    }
  }
}

template < class DomainSpace, class RangeSpace >
template < class GridFunctionType, class DiscreteFunctionType >
void DirichletConstraints<DomainSpace, RangeSpace>::apply( const GridFunctionType& u, DiscreteFunctionType& w ) const
{
  updateDirichletDofs();

  typedef typename DomainSpaceType :: IteratorType IteratorType;
  typedef typename IteratorType :: Entity EntityType;

  // if Dirichlet Dofs have been found, treat them
  if( hasDirichletDofs_ )
  {
    const IteratorType end = domain_space_.end();
    for( IteratorType it = domain_space_.begin(); it != end; ++it )
    {
      const EntityType &entity = *it;
      dirichletDofTreatment( entity, u, w );
    }
  }
}

template < class DomainSpace, class RangeSpace >
template< class LinearOperator, class EntityType >
void DirichletConstraints<DomainSpace, RangeSpace>::dirichletDofsCorrectOnEntity(LinearOperator& linearOperator,
                                                                                 const EntityType &entity ) const
{
  // get slave dof structure (for parallel runs)   /*@LST0S@*/
  SlaveDofsType& slave_dofs = this->slaveDofs();
  const int numSlaveDofs = slave_dofs.size();

  typedef typename DomainSpaceType :: LagrangePointSetType
          LagrangePointSetType;
  const LagrangePointSetType &lagrangePointSet = domain_space_.lagrangePointSet( entity );

  typedef typename LinearOperator :: LocalMatrixType LocalMatrixType;

  // get local matrix from linear operator
  LocalMatrixType localMatrix = linearOperator.localMatrix( entity, entity );

  // get number of basis functions
  const int localBlocks = lagrangePointSet.size();
  const int localBlockSize = DomainSpaceType :: localBlockSize ;

  // map local to global dofs
  std::vector<std::size_t> globalDofs(localBlockSize);
  std::vector<std::size_t> globalBlockDofs(localBlocks);

  domain_space_.mapper().map(entity,globalDofs);
  domain_space_.blockMapper().map(entity,globalBlockDofs);

  // counter for all local dofs (i.e. localBlockDof * localBlockSize + ... )
  int localDof = 0;
  // iterate over face dofs and set unit row
  for( int localBlockDof = 0 ; localBlockDof < localBlocks; ++ localBlockDof )
  {

    if( dirichletBlocks_[ globalBlockDofs[localBlockDof]] )
    {
      for( int l = 0; l < localBlockSize; ++ l, ++ localDof )
      {
        // clear all other columns
        localMatrix.clearRow( localDof );

        // set diagonal to 1
        localMatrix.set( localDof, localDof, 1 );

        // cancel all slave dofs (for parallel runs)
        for( int i = 0; i < numSlaveDofs; ++i )
        {
          // slave dofs are canceled
          if( (int)globalDofs[l] == slave_dofs[ i ] )
            localMatrix.set( localDof, localDof, 0 );
        }
      }
    }
    else
    {
      // increase localDof anyway
      localDof += localBlockSize ;
    }
  }
}

template < class DomainSpace, class RangeSpace >
template< class EntityType, class GridFunctionType, class DiscreteFunctionType >
void DirichletConstraints<DomainSpace, RangeSpace>::dirichletDofTreatment(const EntityType &entity,
                                                                          const GridFunctionType& u,
                                                                          DiscreteFunctionType &w ) const
{
  typedef typename DiscreteFunctionType :: DiscreteFunctionSpaceType DiscreteSpaceType;
  typedef typename GridFunctionType :: LocalFunctionType GridLocalFunctionType;
  typedef typename DiscreteFunctionType :: LocalFunctionType LocalFunctionType;

  typedef typename DiscreteSpaceType :: LagrangePointSetType
                                                     LagrangePointSetType;
  typedef typename DiscreteSpaceType :: GridPartType GridPartType;

  const int faceCodim = 1;
  typedef typename GridPartType :: IntersectionIteratorType
            IntersectionIteratorType;
  typedef typename LagrangePointSetType
  :: template Codim< faceCodim > :: SubEntityIteratorType
            FaceDofIteratorType;

  // get local functions of result
  LocalFunctionType wLocal = w.localFunction( entity );

  // get local functions of argument
  GridLocalFunctionType uLocal = u.localFunction( entity );

  const LagrangePointSetType &lagrangePointSet = domain_space_.lagrangePointSet( entity );

  // get number of Lagrange Points
  const int numBlocks = lagrangePointSet.size();

  int localDof = 0;
  const int localBlockSize = DiscreteSpaceType :: localBlockSize ;

  // map local to global BlockDofs
  std::vector<std::size_t> globalBlockDofs(numBlocks);
  domain_space_.blockMapper().map(entity,globalBlockDofs);

  // iterate over face dofs and set unit row
  for( int localBlock = 0 ; localBlock < numBlocks; ++ localBlock )
  {
    if( dirichletBlocks_[ globalBlockDofs[ localBlock ] ] )
    {
      typedef typename DomainSpaceType :: RangeType RangeType;
      RangeType phi( 0 );

      // evaluate data
      uLocal.evaluate( lagrangePointSet[ localBlock ], phi );

      // store result to dof vector
      for( int l = 0; l < localBlockSize ; ++ l, ++localDof )
      {
        // store result
        wLocal[ localDof ] = phi[ l ];
      }
    }
    else
    {
      // increase localDofs by block size
      localDof += localBlockSize ;
    }

  }
}


template < class DomainSpace, class RangeSpace >
void DirichletConstraints<DomainSpace, RangeSpace>::updateDirichletDofs() const
{
  if( sequence_ != domain_space_.sequence() )
  {
//      // only start search if Dirichlet boundary is present
//      if( ! boundary_.hasDirichletBoundary() )
//      {
//        hasDirichletDofs_ = false ;
//        return ;
//      }

    // resize flag vector with number of blocks and reset flags
    const int blocks = domain_space_.blockMapper().size() ;
    dirichletBlocks_.resize( blocks );
    for( int i=0; i<blocks; ++i )
      dirichletBlocks_[ i ] = false ;

    typedef typename DomainSpaceType :: IteratorType IteratorType;
    typedef typename IteratorType :: Entity EntityType;

    bool hasDirichletBoundary = false;
    const IteratorType end = domain_space_.end();
    for( IteratorType it = domain_space_.begin(); it != end; ++it )
    {
      const EntityType &entity = *it;
      // if entity has boundary intersections
      if( entity.hasBoundaryIntersections() )
      {
        hasDirichletBoundary |= searchEntityDirichletDofs( entity, boundary_ );
      }
    }

    // update sequence number
    sequence_ = domain_space_.sequence();
    hasDirichletDofs_ = hasDirichletBoundary ;
    domain_space_.gridPart().grid().comm().max( hasDirichletDofs_ );
  }
}

template < class DomainSpace, class RangeSpace >
template< class EntityType >
bool DirichletConstraints<DomainSpace, RangeSpace>::searchEntityDirichletDofs(const EntityType &entity,
                                                                              const BoundaryType& boundary) const
{

  typedef typename DomainSpaceType :: LagrangePointSetType
          LagrangePointSetType;

  typedef typename DomainSpaceType :: GridPartType GridPartType;

  const int faceCodim = 1;
  typedef typename GridPartType :: IntersectionIteratorType
            IntersectionIteratorType;

  typedef typename LagrangePointSetType
  :: template Codim< faceCodim > :: SubEntityIteratorType
          FaceDofIteratorType;

  typedef typename DomainSpaceType :: DomainType DomainType ;

  const GridPartType &gridPart = domain_space_.gridPart();

  // default is false
  bool hasDirichletBoundary = false;

  typedef typename EntityType :: Geometry Geometry;
  const Geometry& geo = entity.geometry();

  // get Lagrange pionts from space
  const LagrangePointSetType &lagrangePointSet = domain_space_.lagrangePointSet( entity );

  // get number of Lagrange Points
  const int localBlocks = lagrangePointSet.size();

  //map local to global BlockDofs
  std::vector<size_t> globalBlockDofs(localBlocks);
  // domain_space_.blockMapper().mapEntityDofs(entity,globalBlockDofs);
  domain_space_.blockMapper().map(entity,globalBlockDofs);

  IntersectionIteratorType it = gridPart.ibegin( entity );
  const IntersectionIteratorType endit = gridPart.iend( entity );
  for( ; it != endit; ++it )
  {
    typedef typename IntersectionIteratorType :: Intersection IntersectionType;
    const IntersectionType& intersection = *it;

    // if intersection is with boundary, adjust data
    if( intersection.boundary() )
    {
      // get face number of boundary intersection
      const int face = intersection.indexInInside();

      if (boundary_.dirichlet(intersection))
      {
        // get dof iterators
        FaceDofIteratorType faceIt
                                    = lagrangePointSet.template beginSubEntity< faceCodim >( face );
        const FaceDofIteratorType faceEndIt
                                    = lagrangePointSet.template endSubEntity< faceCodim >( face );
        for( ; faceIt != faceEndIt; ++faceIt )
        {
          // get local dof number (expensive operation, therefore cache result)
          const int localBlock = *faceIt;

          // mark global DoF number
          assert( globalBlockDofs[ localBlock ] < dirichletBlocks_.size() );
          dirichletBlocks_[globalBlockDofs[ localBlock ] ] = true ;

        }
        hasDirichletBoundary = true;
      }
    }
  }

  return hasDirichletBoundary;
}

template < class DomainSpace, class RangeSpace >
typename DirichletConstraints<DomainSpace, RangeSpace>::SlaveDofsType*
DirichletConstraints<DomainSpace, RangeSpace>::getSlaveDofs ( const typename DirichletConstraints<DomainSpace, RangeSpace>::DomainSpaceType &space )
{
  SlaveDofsKeyType key( space, space.mapper() );
  return &(SlaveDofsProviderType :: getObject( key ));
}

template < class DomainSpace, class RangeSpace >
typename DirichletConstraints<DomainSpace, RangeSpace>::SlaveDofsType&
DirichletConstraints<DomainSpace, RangeSpace>::slaveDofs () const
{
  slaveDofs_->rebuild();
  return *slaveDofs_;
}

// explicit instantiations of dirichlet constraints for our function space types
template<> class DirichletConstraints<CommonTraits::DiscreteFunctionSpaceType>;
template<> class DirichletConstraints<MsFEM::MsFEMTraits::SubGridDiscreteFunctionSpaceType>;


void setDirichletDofsZero(CommonTraits::DiscreteFunctionType& function)
{
//  // set dirichlet dofs to zero
//  static const auto boundary = Problem::getModelData()->boundaryInfo();
//  static DirichletConstraints<CommonTraits::DiscreteFunctionSpaceType> constraints(*boundary, function.space());
//  constraints.setValue(0.0, function);
  return;
}

} // namespace Multiscale
} // namespace Dune
