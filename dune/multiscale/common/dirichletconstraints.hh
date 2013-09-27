#ifndef DUNE_DIRICHLETCONSTRAINTS_HH
#define DUNE_DIRICHLETCONSTRAINTS_HH

#include <dune/stuff/grid/boundaryinfo.hh>
#include <dune/stuff/discretefunction/projection/heterogenous.hh>

#include <dune/fem/function/common/scalarproducts.hh>
#include <dune/multiscale/common/traits.hh>
#include <dune/multiscale/msfem/msfem_traits.hh>
#include <dune/multiscale/problems/selector.hh>

namespace Dune {
namespace Multiscale {

template <class DomainSpace, class RangeSpace = DomainSpace>
class DirichletConstraints {
public:
  typedef Dune::Stuff::GridboundaryInterface<typename DomainSpace::GridType::LeafGridView> BoundaryType;
  typedef DomainSpace DomainSpaceType;
  typedef RangeSpace RangeSpaceType;

  //! type of grid partition
  typedef typename DomainSpaceType::GridPartType GridPartType;
  //! type of grid
  typedef typename DomainSpaceType::GridType GridType;

  // types for boundary treatment
  // ----------------------------
  typedef typename DomainSpaceType::MapperType MapperType;
  typedef Fem::SlaveDofs<DomainSpaceType, MapperType> SlaveDofsType;
  typedef typename SlaveDofsType::SingletonKey SlaveDofsKeyType;
  typedef Fem::SingletonList<SlaveDofsKeyType, SlaveDofsType> SlaveDofsProviderType;

  DirichletConstraints(const BoundaryType& boundary, const DomainSpaceType& domain_space)
    : boundary_(boundary)
    , domain_space_(domain_space)
    , slaveDofs_(getSlaveDofs(domain_space_))
    , dirichletBlocks_()
    // mark DoFs on the Dirichlet boundary
    , hasDirichletDofs_(false)
    , sequence_(-1) {}

  /*! treatment of Dirichlet-DoFs for given discrete function
   *
   *   \note A LagrangeDomainSpace is implicitly assumed.
   *
   *   \param[in]  u   discrete function providing the constraints
   *   \param[out] w   discrete function the constraints are applied to
   */
  template <class DiscreteFunctionType>
  void operator()(const DiscreteFunctionType& u, DiscreteFunctionType& w) const {
    updateDirichletDofs();

    // if Dirichlet Dofs have been found, treat them
    if (hasDirichletDofs_) {
      auto uIt = u.dbegin();
      auto wIt = w.dbegin();

      constexpr auto localBlockSize = DiscreteFunctionType::DiscreteFunctionSpaceType::localBlockSize;
      // loop over all blocks
      for (auto blockDof : DSC::valueRange(u.space().blockMapper().size())) {
        if (dirichletBlocks_[blockDof]) {
          // copy dofs of the block
          for (unsigned int l = 0; l < localBlockSize; ++l, ++wIt, ++uIt) {
            assert(uIt != u.dend());
            assert(wIt != w.dend());
            (*wIt) = (*uIt);
          }
        } else {
          // increase dof iterators anyway
          for (unsigned int l = 0; l < localBlockSize; ++l, ++wIt, ++uIt) {
          }
        }
      }
    }
  }

  /*! treatment of Dirichlet-DoFs for given discrete function
 *
 *   \note A LagrangeDomainSpace is implicitly assumed.
 *
 *   \param[in]  val a value that will be used to set all dirichlet dofs
 *   \param[out] w   discrete function the constraints are applied to
 */
  template <class DiscreteFunctionType>
  void setValue(const typename DiscreteFunctionType::RangeFieldType val, DiscreteFunctionType& w) const {
    updateDirichletDofs();

    // if Dirichlet Dofs have been found, treat them
    if (hasDirichletDofs_) {
      auto wIt = w.dbegin();

      const auto localBlockSize = DiscreteFunctionType::DiscreteFunctionSpaceType::localBlockSize;
      // loop over all blocks
      const auto blocks = w.space().blockMapper().size();
      for (unsigned int blockDof = 0; blockDof < blocks; ++blockDof) {
        if (dirichletBlocks_[blockDof]) {
          // copy dofs of the block
          for (unsigned int l = 0; l < localBlockSize; ++l, ++wIt) {
            assert(wIt != w.dend());
            (*wIt) = val;
          }
        } else {
          // increase dof iterators anyway
          for (unsigned int l = 0; l < localBlockSize; ++l, ++wIt) {
          }
        }
      }
    }
  }

  /*! treatment of Dirichlet-DoFs for given discrete function
   *
   *   \note A LagrangeDomainSpace is implicitly assumed.
   *
   *   \param[in]  u   discrete function providing the constraints
   *   \param[out] w   discrete function the constraints are applied to
   */
  template <class GridFunctionType, class DiscreteFunctionType>
  void operator()(const GridFunctionType& u, DiscreteFunctionType& w) const {
    apply(u, w);
  }

  /*! treatment of Dirichlet-DoFs for solution and right-hand-side
   *
   *   delete rows for dirichlet-DoFs, setting diagonal element to 1.
   *
   *   \note A LagrangeDomainSpace is implicitly assumed.
   *
   *   \param[out] linearOperator  linear operator to be adjusted
   */
  template <class LinearOperator>
  void applyToOperator(LinearOperator& linearOperator) const {
    updateDirichletDofs();
    // if Dirichlet Dofs have been found, treat them
    if (hasDirichletDofs_) {
      for (const auto& entity : domain_space_) {
        // adjust linear operator
        dirichletDofsCorrectOnEntity(linearOperator, entity);
      }
    }
  }

protected:
  template <class GridFunctionType, class DiscreteFunctionType>
  void apply(const GridFunctionType& u, DiscreteFunctionType& w) const {
    updateDirichletDofs();
    // if Dirichlet Dofs have been found, treat them
    if (hasDirichletDofs_) {
      for (const auto& entity : domain_space_) {
        dirichletDofTreatment(entity, u, w);
      }
    }
  }

  /*! treatment of Dirichlet-DoFs for one entity
   *
   *   delete rows for dirichlet-DoFs, setting diagonal element to 1.
   *
   *   \note A LagrangeDomainSpace is implicitly assumed.
   *
   *   \param[in]  entity  entity to perform Dirichlet treatment on
   */
  template <class LinearOperator, class EntityType>
  void dirichletDofsCorrectOnEntity(LinearOperator& linearOperator, const EntityType& entity) const {
    // get slave dof structure (for parallel runs)   /*@LST0S@*/
    auto& slave_dofs = this->slaveDofs();
    const auto numSlaveDofs = slave_dofs.size();
    const auto& lagrangePointSet = domain_space_.lagrangePointSet(entity);
    auto localMatrix = linearOperator.localMatrix(entity, entity);

    // get number of basis functions
    const auto localBlocks = lagrangePointSet.size();
    const auto localBlockSize = DomainSpaceType::localBlockSize;

    // map local to global dofs
    std::vector<std::size_t> globalDofs(localBlockSize);
    std::vector<std::size_t> globalBlockDofs(localBlocks);

    domain_space_.mapper().map(entity, globalDofs);
    domain_space_.blockMapper().map(entity, globalBlockDofs);

    // counter for all local dofs (i.e. localBlockDof * localBlockSize + ... )
    int localDof = 0;
    // iterate over face dofs and set unit row
    for (auto localBlockDof : DSC::valueRange(localBlocks)) {

      if (dirichletBlocks_[globalBlockDofs[localBlockDof]]) {
        for (int l = 0; l < localBlockSize; ++l, ++localDof) {
          // clear all other columns
          localMatrix.clearRow(localDof);

          // set diagonal to 1
          localMatrix.set(localDof, localDof, 1);

          // cancel all slave dofs (for parallel runs)
          for (int i = 0; i < numSlaveDofs; ++i) {
            // slave dofs are canceled
            if ((int)globalDofs[l] == slave_dofs[i])
              localMatrix.set(localDof, localDof, 0);
          }
        }
      } else {
        // increase localDof anyway
        localDof += localBlockSize;
      }
    }
  }

  //! set the dirichlet points to exact values
  template <class EntityType, class GridFunctionType, class DiscreteFunctionType>
  void dirichletDofTreatment(const EntityType& entity, const GridFunctionType& u, DiscreteFunctionType& w) const {
    typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteSpaceType;

    auto wLocal = w.localFunction(entity);
    auto uLocal = u.localFunction(entity);
    const auto& lagrangePointSet = domain_space_.lagrangePointSet(entity);
    const auto numBlocks = lagrangePointSet.size();

    int localDof = 0;
    const auto localBlockSize = DiscreteSpaceType::localBlockSize;

    // map local to global BlockDofs
    std::vector<std::size_t> globalBlockDofs(numBlocks);
    domain_space_.blockMapper().map(entity, globalBlockDofs);

    // iterate over face dofs and set unit row
    for (auto localBlock : DSC::valueRange(numBlocks)) {
      if (dirichletBlocks_[globalBlockDofs[localBlock]]) {
        typedef typename DomainSpaceType::RangeType RangeType;
        RangeType phi(0);

        // evaluate data
        uLocal.evaluate(lagrangePointSet[localBlock], phi);

        // store result to dof vector
        for (int l = 0; l < localBlockSize; ++l, ++localDof) {
          // store result
          wLocal[localDof] = phi[l];
        }
      } else {
        // increase localDofs by block size
        localDof += localBlockSize;
      }
    }
  }

protected:
  // detect all DoFs on the Dirichlet boundary
  void updateDirichletDofs() const {
    if (sequence_ != domain_space_.sequence()) {
      //      // only start search if Dirichlet boundary is present
      //      if( ! boundary_.hasDirichletBoundary() )
      //      {
      //        hasDirichletDofs_ = false ;
      //        return ;
      //      }

      // resize flag vector with number of blocks and reset flags
      const auto blocks = domain_space_.blockMapper().size();
      dirichletBlocks_ = std::vector<bool>(blocks, false);

      bool hasDirichletBoundary = false;
      for (const auto& entity : domain_space_) {
        // if entity has boundary intersections
        if (entity.hasBoundaryIntersections()) {
          hasDirichletBoundary |= searchEntityDirichletDofs(entity, boundary_);
        }
      }

      // update sequence number
      sequence_ = domain_space_.sequence();
      hasDirichletDofs_ = hasDirichletBoundary;
      domain_space_.gridPart().grid().comm().max(hasDirichletDofs_);
    }
  }

  // detect all DoFs on the Dirichlet boundary of the given entity
  template <class EntityType>
  bool searchEntityDirichletDofs(const EntityType& entity, const BoundaryType& /*boundary*/) const {
    static const int faceCodim = 1;
    const auto& gridPart = domain_space_.gridPart();
    bool hasDirichletBoundary = false;
    const auto& lagrangePointSet = domain_space_.lagrangePointSet(entity);
    const auto localBlocks = lagrangePointSet.size();

    // map local to global BlockDofs
    std::vector<size_t> globalBlockDofs(localBlocks);
    // domain_space_.blockMapper().mapEntityDofs(entity,globalBlockDofs);
    domain_space_.blockMapper().map(entity, globalBlockDofs);

    for (const auto& intersection : DSC::intersectionRange(gridPart, entity)) {
      // if intersection is with boundary, adjust data
      if (intersection.boundary()) {
        // get face number of boundary intersection
        const int face = intersection.indexInInside();

        if (boundary_.dirichlet(intersection)) {
          // get dof iterators
          for (const auto& lp : DSC::lagrangePointSetRange(lagrangePointSet, face)) {
            // mark global DoF number
            assert(globalBlockDofs[lp] < dirichletBlocks_.size());
            dirichletBlocks_[globalBlockDofs[lp]] = true;
          }
          hasDirichletBoundary = true;
        }
      }
    }

    return hasDirichletBoundary;
  }

  //! pointer to slave dofs
  const BoundaryType& boundary_;
  const DomainSpaceType& domain_space_;
  SlaveDofsType* const slaveDofs_;
  mutable std::vector<bool> dirichletBlocks_;
  mutable bool hasDirichletDofs_;
  mutable int sequence_;

  // return slave dofs
  static SlaveDofsType* getSlaveDofs(const DomainSpaceType& space) {
    SlaveDofsKeyType key(space, space.mapper());
    return &(SlaveDofsProviderType::getObject(key));
  }

  // return reference to slave dofs
  SlaveDofsType& slaveDofs() const {
    slaveDofs_->rebuild();
    return *slaveDofs_;
  }
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
DirichletConstraints<CommonTraits::DiscreteFunctionSpaceType>&
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
DirichletConstraints<CommonTraits::DiscreteFunctionSpaceType>&
getConstraintsFine(const CommonTraits::DiscreteFunctionSpaceType& space);

void projectDirichletValues(const CommonTraits::DiscreteFunctionSpaceType& coarseSpace,
                            CommonTraits::DiscreteFunctionType& function);

} // end namespace Multiscale
} // end namespace Dune
#endif
