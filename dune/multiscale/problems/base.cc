#include <config.h>
#include "base.hh"

void Dune::Multiscale::Problem::DiffusionBase::jacobianDiffusiveFlux(const Dune::Multiscale::Problem::DomainType& /*x*/, const Dune::Multiscale::Problem::JacobianRangeType &, const Dune::Multiscale::Problem::JacobianRangeType& /*direction_gradient*/, Dune::Multiscale::Problem::JacobianRangeType& /*flux*/) const
{
  DUNE_THROW(NotImplemented, "");
}
