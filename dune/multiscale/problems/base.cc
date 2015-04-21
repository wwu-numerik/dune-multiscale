#include <config.h>
#include "base.hh"

void DMP::DiffusionBase::jacobianDiffusiveFlux(
    const DMP::DomainType& /*x*/, const DMP::JacobianRangeType&,
    const DMP::JacobianRangeType& /*direction_gradient*/,
    DMP::JacobianRangeType& /*flux*/) const {
  DUNE_THROW(NotImplemented, "");
}
