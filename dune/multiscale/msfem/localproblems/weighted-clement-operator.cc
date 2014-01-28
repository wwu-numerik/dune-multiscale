#include <config.h>
#include <assert.h>
#include <dune/common/timer.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/stuff/fem/localmatrix_proxy.hh>
#include <dune/stuff/common/ranges.hh>
#include <dune/stuff/aliases.hh>
#include <dune/multiscale/tools/misc/linear-lagrange-interpolation.hh>
#include <dune/stuff/fem/functions/integrals.hh>
#include "weighted-clement-operator.hh"

namespace Dune {
namespace Multiscale {
namespace MsFEM {
class MacroMicroGridSpecifier;
}  // namespace MsFEM
}  // namespace Multiscale
}  // namespace Dune

DMM::WeightedClementOperator::WeightedClementOperator(
    const DMM::WeightedClementOperator::DiscreteFunctionSpaceType& space,
    const DMM::WeightedClementOperator::CoarseDiscreteFunctionSpaceType& coarse_space,
    const DMM::WeightedClementOperator::CoarseNodeVectorType& coarse_nodes,
    const DMM::WeightedClementOperator::CoarseBasisFunctionList& coarse_basis,
    const std::map<std::size_t, std::size_t>& global_id_to_internal_id,
    const DMM::MacroMicroGridSpecifier& specifier)
  : discreteFunctionSpace_(space)
  , coarse_space_(coarse_space)
  , dofManager_(DofManagerType::instance(space.grid()))
  , specifier_(specifier)
  , sparsity_pattern_(discreteFunctionSpace_, coarse_space_, specifier_)
  , linearOperator_(discreteFunctionSpace_, coarse_space_)
  , coarse_nodes_(coarse_nodes)
  , coarse_basis_(coarse_basis)
  , global_id_to_internal_id_(global_id_to_internal_id)
  , sequence_(-1)
  , gradCache_(discreteFunctionSpace_.mapper().maxNumDofs())
  , values_(discreteFunctionSpace_.mapper().maxNumDofs()) {}

void DMM::WeightedClementOperator::
operator()(const DMM::WeightedClementOperator::DiscreteFunctionType& u,
           DMM::WeightedClementOperator::CoarseDiscreteFunctionType& w) const {
  systemMatrix().apply(u, w); /*@\label{sto:matrixEval}@*/
}

const DMM::WeightedClementOperator::PreconditionMatrixType&
DMM::WeightedClementOperator::preconditionMatrix() const {
  return systemMatrix().preconditionMatrix();
}

void DMM::WeightedClementOperator::applyTransposed(
    const DMM::WeightedClementOperator::CoarseDiscreteFunctionType& u,
    DMM::WeightedClementOperator::DiscreteFunctionType& w) const {
  systemMatrix().apply_t(u, w); /*@\label{sto:applytransposed}@*/
}

bool DMM::WeightedClementOperator::hasPreconditionMatrix() const {
  return linearOperator_.hasPreconditionMatrix();
}

void DMM::WeightedClementOperator::print(std::ostream& out) const {
  systemMatrix().matrix().print(out);
}

const DMM::WeightedClementOperator::DiscreteFunctionSpaceType&
DMM::WeightedClementOperator::discreteFunctionSpace() const {
  return discreteFunctionSpace_;
}

const DMM::WeightedClementOperator::LinearOperatorType&
DMM::WeightedClementOperator::systemMatrix() const {
  // if stored sequence number it not equal to the one of the
  // dofManager (or space) then the grid has been changed
  // and matrix has to be assembled new
  if (sequence_ != dofManager_.sequence()) /*@\label{sto:sequence}@*/
    assemble();

  return linearOperator_;
}

void DMM::WeightedClementOperator::boundaryTreatment() const {
  for (const auto& entity : discreteFunctionSpace_) {
    for (const auto& coarse_entity : coarse_space_) {
      if (!DSG::entities_identical(entity, coarse_entity))
        continue;

      // if entity has boundary intersections
      if (entity.hasBoundaryIntersections()) {
        // get local matrix from matrix object
        DSFe::LocalMatrixProxy<LinearOperatorType> localMatrix(linearOperator_, entity, coarse_entity);

        const auto& lagrangePointSet = discreteFunctionSpace_.lagrangePointSet(entity);

        const auto endiit = discreteFunctionSpace_.gridPart().iend(entity);
        for (auto iit = discreteFunctionSpace_.gridPart().ibegin(entity); iit != endiit; ++iit) {
          if (iit->neighbor()) // if there is a neighbor entity
            continue;

          const int face = (*iit).indexInInside();
          for (const auto& lp : DSC::lagrangePointSetRange(lagrangePointSet, face))
            localMatrix.unitRow(lp);
        }
      }
    }
  }
}
