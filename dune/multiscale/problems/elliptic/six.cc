#include <config.h>
#include <assert.h>
#include <dune/common/exceptions.hh>
#include <dune/stuff/common/parameter/validation.hh>
#include <sstream>

#include "dune/multiscale/problems/base.hh"
#include "dune/multiscale/problems/constants.hh"
#include "six.hh"

namespace Dune {
namespace Multiscale {
namespace Problem {
namespace Six {

// default value for epsilon (if not sprecified in the parameter file)
CONSTANTSFUNCTION(0.05)

ModelProblemData::ModelProblemData() : IModelProblemData(constants()) {
  assert(constants_.epsilon != 0.0);
  if (constants().get("stochastic_pertubation", false) && !(this->problemAllowsStochastics()))
    DUNE_THROW(Dune::InvalidStateException,
               "The problem does not allow stochastic perturbations. Please, switch the key off.");
}

std::string ModelProblemData::getMacroGridFile() const {
  return ("../dune/multiscale/grids/macro_grids/elliptic/cube_two.dgf");
}

bool ModelProblemData::problemIsPeriodic() const {
  return true; // = problem is periodic
}

bool ModelProblemData::problemAllowsStochastics() const {
  return true; // = problem allows stochastic perturbations
}

void Diffusion::diffusiveFlux(
    const DomainType& x,
    const JacobianRangeType& gradient,
    JacobianRangeType& flux) const {

  // coeff.first = 1.01 + cos( 2.0 * M_PI * (x[0] / epsilon) ) + stochastic perturbation
  // coeff.second = 1.01 + cos( 2.0 * M_PI * (x[0] / epsilon) ) + stochastic perturbation
  const auto coeff = constants().coefficients(x);
  flux[0][0] = coeff.first * gradient[0][0];
  flux[0][1] = coeff.second * gradient[0][1];
}

void Diffusion::jacobianDiffusiveFlux(
    const DomainType& x,
    const JacobianRangeType&,
    const JacobianRangeType& direction_gradient,
    JacobianRangeType& flux) const {

  // coeff.first = 1.01 + cos( 2.0 * M_PI * (x[0] / epsilon) ) + stochastic perturbation
  // coeff.second = 1.01 + cos( 2.0 * M_PI * (x[0] / epsilon) ) + stochastic perturbation
  const auto coeff = constants().coefficients(x);

  flux[0][0] = coeff.first * direction_gradient[0][0];
  flux[0][1] = coeff.second * direction_gradient[0][1];
}

void ExactSolution::evaluate(
    const ExactSolution::DomainType&,
    ExactSolution::RangeType&) const {
  DUNE_THROW(Dune::NotImplemented, "Exact solution not available!");
}

void ExactSolution::jacobian(
    const ExactSolution::DomainType&,
    ExactSolution::JacobianRangeType&) const {
  DUNE_THROW(Dune::NotImplemented, "Exact solution not available!");
}

void ExactSolution::evaluate(
    const ExactSolution::DomainType& x, const TimeType&,
    ExactSolution::RangeType& y) const {
  evaluate(x, y);
}

} // namespace Six
} // namespace Problem
} // namespace Multiscale {
} // namespace Dune {
