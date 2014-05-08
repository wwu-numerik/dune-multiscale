#ifndef DUNE_DIRICHLETCONSTRAINTS_HH
#define DUNE_DIRICHLETCONSTRAINTS_HH

#include <assert.h>
#include <dune/fem/function/common/scalarproducts.hh>
#include <dune/fem/storage/singletonlist.hh>
#include <dune/multiscale/common/traits.hh>
#include <dune/multiscale/msfem/msfem_traits.hh>
#include <dune/multiscale/problems/selector.hh>
#include <dune/stuff/common/ranges.hh>
#include <dune/stuff/discretefunction/projection/heterogenous.hh>
#include <dune/stuff/grid/boundaryinfo.hh>
#include <dune/stuff/functions/femadapter.hh>
#include <cstddef>
#include <vector>

namespace Dune {
namespace Multiscale {

/** explicitly instantiated in dirichletconstraints.cc for Coarse/Local DF types
 **/
template <class DiscreteFunctionImp>
class DirichletConstraints {

  typedef DiscreteFunctionImp DiscreteFunctionType;
  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DomainSpaceType;
  typedef DSG::BoundaryInfoInterface<typename DomainSpaceType::GridType::LeafGridView::Intersection>
  BoundaryType;
  typedef typename DomainSpaceType::GridPartType GridPartType;
  typedef typename DomainSpaceType::GridType GridType;

  typedef typename BackendChooser<DomainSpaceType>::LinearOperatorType LinearOperatorType;
  typedef typename DomainSpaceType::EntityType EntityType;

  typedef Dune::Fem::GridFunctionAdapter<CommonTraits::DirichletDataType,
                                         typename DiscreteFunctionImp::GridPartType> GridFunctionType;

public:
  DirichletConstraints(const BoundaryType& boundary, const DomainSpaceType& domain_space);

private:
  /*! treatment of Dirichlet-DoFs for given discrete function
   *
   *   \note A LagrangeDomainSpace is implicitly assumed.
   *
   *   \param[in]  u   discrete function providing the constraints
   *   \param[out] w   discrete function the constraints are applied to
   */
  void operator()(const DiscreteFunctionType& u, DiscreteFunctionType& w) const;

public:
  /*! treatment of Dirichlet-DoFs for given discrete function
   *
   *   \note A LagrangeDomainSpace is implicitly assumed.
   *
   *   \param[in]  u   discrete function providing the constraints
   *   \param[out] w   discrete function the constraints are applied to
   */
  void operator()(const GridFunctionType& u, DiscreteFunctionType& w) const;

  /*! treatment of Dirichlet-DoFs for given discrete function
 *
 *   \note A LagrangeDomainSpace is implicitly assumed.
 *
 *   \param[in]  val a value that will be used to set all dirichlet dofs
 *   \param[out] w   discrete function the constraints are applied to
 */
  void setValue(const typename DiscreteFunctionType::RangeFieldType val, DiscreteFunctionType& w) const;

  /*! treatment of Dirichlet-DoFs for solution and right-hand-side
   *
   *   delete rows for dirichlet-DoFs, setting diagonal element to 1.
   *
   *   \note A LagrangeDomainSpace is implicitly assumed.
   *
   *   \param[out] linearOperator  linear operator to be adjusted
   */
  void applyToOperator(LinearOperatorType& linearOperator) const;

protected:
  void apply(const GridFunctionType& u, DiscreteFunctionType& w) const;

  /*! treatment of Dirichlet-DoFs for one entity
   *
   *   delete rows for dirichlet-DoFs, setting diagonal element to 1.
   *
   *   \note A LagrangeDomainSpace is implicitly assumed.
   *
   *   \param[in]  entity  entity to perform Dirichlet treatment on
   */
  void dirichletDofsCorrectOnEntity(LinearOperatorType& linearOperator, const EntityType& entity) const;

  //! set the dirichlet points to exact values
  void dirichletDofTreatment(const EntityType& entity, const GridFunctionType& u, DiscreteFunctionType& w) const;

  //! detect all DoFs on the Dirichlet boundary
  void updateDirichletDofs() const;

  //! detect all DoFs on the Dirichlet boundary of the given entity
  bool searchEntityDirichletDofs(const EntityType& entity, const BoundaryType& /*boundary*/) const;

  const BoundaryType& boundary_;
  const DomainSpaceType& domain_space_;
  mutable std::vector<bool> dirichletBlocks_;
  mutable bool hasDirichletDofs_;
  mutable int sequence_;
};

/** Get the constraints for a given discrete function space.
*
* @param space The discrete function space.
*
* @attention As the dirichlet constraints are stored in a static variable in
*this function, probably something
*            bad will happen if you run this function with spaces of the same
*type living on different grid(part)
 *           instances. Therefore, this method should only be used for spaces on
*the coarse grid!
*/
DirichletConstraints<CommonTraits::DiscreteFunctionType>&
getConstraintsCoarse(const CommonTraits::DiscreteFunctionSpaceType& space);

/** Get the constraints for a given discrete function space.
*
* @param space The discrete function space.
*
* @attention As the dirichlet constraints are stored in a static variable in
*this function, probably something
*            bad will happen if you run this function with spaces of the same
*type living on different grid(part)
 *           instances. Therefore, this method should only be used for spaces on
*the fine grid!
*/
DirichletConstraints<CommonTraits::DiscreteFunctionType>&
getConstraintsFine(const CommonTraits::DiscreteFunctionSpaceType& space);

/** Copy the dirichlet values to a given discrete function on the fine mesh.
 *
 *  This method projects the dirichlet values to a function on the coarse mesh.
 *  In a second step, that function is projected to a given function on the fine mesh.
 *
 * @attention The given function will be cleared!
 * @note suboptimal since local_grid change
 *
 * @param[in] coarseSpace the discrete function space on the coarse grid, needed
 *                        to perfom the correct projection.
 * @param[in, out] function The function to which the values will be projected.
 */
template <class LocalGridOrCoarseDiscreteFunctionSpaceImp>
void copyDirichletValues(const CommonTraits::DiscreteFunctionSpaceType& coarseSpace,
                         Fem::DiscreteFunctionInterface<LocalGridOrCoarseDiscreteFunctionSpaceImp>& function) {
  static bool initialized = false;
  static CommonTraits::DiscreteFunctionType dirichletExtensionCoarse("Dirichlet Extension Coarse", coarseSpace);
  if (!initialized) {
    initialized = true;

    const auto& dirichletDataPtr = Problem::getDirichletData();
    const auto& dirichletData = *dirichletDataPtr;

    const auto gf = DS::gridFunctionAdapter(dirichletData, coarseSpace.gridPart());
    dirichletExtensionCoarse.clear();
    getConstraintsCoarse(coarseSpace)(gf, dirichletExtensionCoarse);
  }
  Dune::Stuff::HeterogenousProjection<>::project(dirichletExtensionCoarse, function);
}

} // end namespace Multiscale
} // end namespace Dune

#endif
