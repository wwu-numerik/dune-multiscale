#include <config.h>
#include <dune/fem/function/common/gridfunctionadapter.hh>
#include <dune/stuff/discretefunction/projection/heterogenous.hh>
#include <memory>

#include "dirichletconstraints.hh"
#include "dune/multiscale/common/traits.hh"
#include "dune/multiscale/problems/base.hh"
#include "dune/multiscale/problems/selector.hh"
#include <dune/multiscale/msfem/msfem_traits.hh>

namespace Dune {
namespace Multiscale {

DirichletConstraints<CommonTraits::DiscreteFunctionType>&
getConstraintsCoarse(const CommonTraits::DiscreteFunctionSpaceType& space) {
  // set dirichlet dofs to zero
  static const auto& boundary = Problem::getModelData()->boundaryInfo();
  static DirichletConstraints<CommonTraits::DiscreteFunctionType> constraints(boundary, space);
  return constraints;
}

DirichletConstraints<CommonTraits::DiscreteFunctionType>&
getConstraintsFine(const CommonTraits::DiscreteFunctionSpaceType& space) {
  // set dirichlet dofs to zero
  static const auto& boundary = Problem::getModelData()->boundaryInfo();
  static DirichletConstraints<CommonTraits::DiscreteFunctionType> constraints(boundary, space);
  return constraints;
}

template <class DF>
DirichletConstraints<DF>::DirichletConstraints(const DirichletConstraints::BoundaryType& boundary,
                                               const DirichletConstraints::DomainSpaceType& domain_space)
  : boundary_(boundary)
  , domain_space_(domain_space)
  , dirichletBlocks_()
  , hasDirichletDofs_(false)
  , sequence_(-1) {}

template <class DF>
void DirichletConstraints<DF>::operator()(const DirichletConstraints::GridFunctionType& u,
                                          DirichletConstraints::DiscreteFunctionType& w) const {
  apply(u, w);
}

template <class DF>
void DirichletConstraints<DF>::setValue(const typename DiscreteFunctionType::RangeFieldType val,
                                        DiscreteFunctionType& w) const {
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

template <class DF>
void DirichletConstraints<DF>::applyToOperator(DirichletConstraints::LinearOperatorType& linearOperator) const {
  updateDirichletDofs();
  // if Dirichlet Dofs have been found, treat them
  if (hasDirichletDofs_) {
    for (const auto& entity : domain_space_) {
      // adjust linear operator
      dirichletDofsCorrectOnEntity(linearOperator, entity);
    }
  }
  linearOperator.communicate();
}

template <class DF>
void DirichletConstraints<DF>::apply(const DirichletConstraints::GridFunctionType& u, DiscreteFunctionType& w) const {
  updateDirichletDofs();
  // if Dirichlet Dofs have been found, treat them
  if (hasDirichletDofs_) {
    for (const auto& entity : domain_space_) {
      dirichletDofTreatment(entity, u, w);
    }
  }
}

template <class DF>
void DirichletConstraints<DF>::dirichletDofsCorrectOnEntity(DirichletConstraints::LinearOperatorType& linearOperator,
                                                            const EntityType& entity) const {
  const auto& lagrangePointSet = domain_space_.lagrangePointSet(entity);
  auto localMatrix = linearOperator.localMatrix(entity, entity);

  // get number of basis functions
  const auto localBlocks = lagrangePointSet.size();
  const auto localBlockSize = DomainSpaceType::localBlockSize;

  // map local to global dofs
  std::vector<std::size_t> globalDofs(localBlockSize);
  std::vector<std::size_t> globalBlockDofs(localBlocks);

  domain_space_.blockMapper().map(entity, globalDofs);
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
      }
    } else {
      // increase localDof anyway
      localDof += localBlockSize;
    }
  }
}

template <class DF>
void DirichletConstraints<DF>::dirichletDofTreatment(const EntityType& entity,
                                                     const DirichletConstraints::GridFunctionType& u,
                                                     DiscreteFunctionType& w) const {
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

template <class DF>
void DirichletConstraints<DF>::updateDirichletDofs() const {
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
//    domain_space_.gridPart().grid().comm().max(hasDirichletDofs_);
  }
}

template <class DF>
bool DirichletConstraints<DF>::searchEntityDirichletDofs(const EntityType& entity,
                                                         const DirichletConstraints::BoundaryType&) const {
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

template <class DF>
void DirichletConstraints<DF>::operator()(const DiscreteFunctionType& u, DiscreteFunctionType& w) const {
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

template class DirichletConstraints<CommonTraits::DiscreteFunctionType>;
//template class DirichletConstraints<DMM::MsFEMTraits::LocalGridDiscreteFunctionType>;

} // namespace Multiscale
} // namespace Dune
