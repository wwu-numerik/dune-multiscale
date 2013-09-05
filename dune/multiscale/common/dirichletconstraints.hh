#ifndef DUNE_DIRICHLETCONSTRAINTS_HH
#define DUNE_DIRICHLETCONSTRAINTS_HH

#include <dune/stuff/grid/boundaryinfo.hh>
#include <dune/fem/function/common/scalarproducts.hh>

namespace Dune { 
namespace Multiscale {

template < class DomainSpace, class RangeSpace = DomainSpace >
class DirichletConstraints 
{
public:
  typedef Dune::Stuff::GridboundaryInterface<typename DomainSpace::GridType::LeafGridView> BoundaryType;
  typedef DomainSpace DomainSpaceType;
  typedef RangeSpace RangeSpaceType;

  //! type of grid partition
  typedef typename DomainSpaceType :: GridPartType GridPartType;
  //! type of grid
  typedef typename DomainSpaceType :: GridType GridType;

  // types for boundary treatment
  // ----------------------------
  typedef typename DomainSpaceType :: MapperType MapperType;
  typedef Fem::SlaveDofs< DomainSpaceType, MapperType > SlaveDofsType;
  typedef typename SlaveDofsType :: SingletonKey SlaveDofsKeyType; 
  typedef Fem::SingletonList< SlaveDofsKeyType, SlaveDofsType > SlaveDofsProviderType;

  DirichletConstraints( const BoundaryType &boundary, const DomainSpaceType& domain_space );

  /*! treatment of Dirichlet-DoFs for given discrete function 
   *
   *   \note A LagrangeDomainSpace is implicitly assumed.
   *
   *   \param[in]  u   discrete function providing the constraints 
   *   \param[out] w   discrete function the constraints are applied to
   */
  template < class DiscreteFunctionType >
  void operator ()( const DiscreteFunctionType& u, DiscreteFunctionType& w ) const;

  /*! treatment of Dirichlet-DoFs for given discrete function
 *
 *   \note A LagrangeDomainSpace is implicitly assumed.
 *
 *   \param[in]  val a value that will be used to set all dirichlet dofs
 *   \param[out] w   discrete function the constraints are applied to
 */
  template < class DiscreteFunctionType >
  void setValue(const typename DiscreteFunctionType::RangeFieldType val, DiscreteFunctionType& w ) const;

  /*! treatment of Dirichlet-DoFs for given discrete function 
   *
   *   \note A LagrangeDomainSpace is implicitly assumed.
   *
   *   \param[in]  u   discrete function providing the constraints 
   *   \param[out] w   discrete function the constraints are applied to
   */
  template < class GridFunctionType, class DiscreteFunctionType >
  void operator ()( const GridFunctionType& u, DiscreteFunctionType& w ) const;

  /*! treatment of Dirichlet-DoFs for solution and right-hand-side
   *
   *   delete rows for dirichlet-DoFs, setting diagonal element to 1.
   *
   *   \note A LagrangeDomainSpace is implicitly assumed.
   *
   *   \param[out] linearOperator  linear operator to be adjusted 
   */
  template <class LinearOperator>
  void applyToOperator( LinearOperator& linearOperator ) const;

protected:  
  template < class GridFunctionType, class DiscreteFunctionType >
  void apply( const GridFunctionType& u, DiscreteFunctionType& w ) const;
  /*! treatment of Dirichlet-DoFs for one entity
   *
   *   delete rows for dirichlet-DoFs, setting diagonal element to 1.
   *
   *   \note A LagrangeDomainSpace is implicitly assumed.
   *
   *   \param[in]  entity  entity to perform Dirichlet treatment on
   */
  template< class LinearOperator, class EntityType >                           
  void dirichletDofsCorrectOnEntity ( LinearOperator& linearOperator, 
                                      const EntityType &entity ) const;

  //! set the dirichlet points to exact values
  template< class EntityType, class GridFunctionType, class DiscreteFunctionType >
  void dirichletDofTreatment( const EntityType &entity,
                              const GridFunctionType& u, 
                              DiscreteFunctionType &w ) const;

protected:
  // detect all DoFs on the Dirichlet boundary 
  void updateDirichletDofs() const;

  // detect all DoFs on the Dirichlet boundary of the given entity 
  template< class EntityType > 
  bool searchEntityDirichletDofs( const EntityType &entity, const BoundaryType& boundary ) const;

  // return slave dofs
  static SlaveDofsType *getSlaveDofs ( const DomainSpaceType &space );

  // return reference to slave dofs
  SlaveDofsType &slaveDofs () const;

  //! pointer to slave dofs 
  const BoundaryType& boundary_;
  const DomainSpaceType& domain_space_;
  SlaveDofsType *const slaveDofs_;
  mutable std::vector< bool > dirichletBlocks_;
  mutable bool hasDirichletDofs_ ;
  mutable int sequence_ ;
};


} // end namespace Multiscale
} // end namespace Dune
#endif
