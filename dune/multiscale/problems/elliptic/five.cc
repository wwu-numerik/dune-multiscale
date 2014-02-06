#include <config.h>
#include <assert.h>
#include <dune/common/exceptions.hh>
#include <dune/stuff/common/parameter/validation.hh>
#include <sstream>

#include "dune/multiscale/problems/base.hh"
#include "dune/multiscale/problems/constants.hh"
#include "five.hh"

namespace Dune {
namespace Multiscale {
namespace Problem {
namespace Five {

// default value for epsilon (if not sprecified in the parameter file)
CONSTANTSFUNCTION(0.05)

ModelProblemData::ModelProblemData() : IModelProblemData(constants()) {
  assert(constants_.epsilon != 0.0);
  if (constants().get("stochastic_pertubation", false) && !(this->problemAllowsStochastics()))
    DUNE_THROW(Dune::InvalidStateException,
               "The problem does not allow stochastic perturbations. Please, switch the key off.");
}

std::string ModelProblemData::getMacroGridFile() const {
  return ("../dune/multiscale/grids/macro_grids/elliptic/corner_singularity.dgf");
}

bool ModelProblemData::problemIsPeriodic() const {
  return true; // = problem is periodic
}

bool ModelProblemData::problemAllowsStochastics() const {
  return true; // = problem allows stochastic perturbations
}

void FirstSource::evaluate(const FirstSource::DomainType& x, FirstSource::RangeType& y) const {

  // circle of radius 0.2 around the reentrant corner at (0.5,0.5)
  double distance = sqrt(pow(x[0] - 0.5, 2.0) + pow(x[1] - 0.5, 2.0));

  if (distance < 0.2) {
    y = 1.0;
  } else {
    y = 0.1;
  }
}

void FirstSource::evaluate(const FirstSource::DomainType& x, const TimeType&, FirstSource::RangeType& y) const {
  evaluate(x, y);
}

void Diffusion::diffusiveFlux(const DomainType& x, const JacobianRangeType& gradient, JacobianRangeType& flux) const {

  // coeff.first = ( 0.1 + ( 1.0 * pow(cos( 2.0 * M_PI * (x[0] / epsilon) ), 2.0) ) ) + stochastic perturbation
  // coeff.second = ( 0.1 + 1e-3 + ( 0.1 * sin( 2.0 * M_PI * (x[1] / epsilon) ) ) ) + stochastic perturbation
  const auto coeff = constants().coefficients_variant_A(x);

  flux[0][0] = coeff.first * (gradient[0][0] + ((1.0 / 3.0) * pow(gradient[0][0], 3.0)));
  flux[0][1] = coeff.second * (gradient[0][1] + ((1.0 / 3.0) * pow(gradient[0][1], 3.0)));
}

void Diffusion::jacobianDiffusiveFlux(const DomainType& x, const JacobianRangeType& position_gradient,
                                      const JacobianRangeType& direction_gradient, JacobianRangeType& flux) const {

  // coeff.first = ( 0.1 + ( 1.0 * pow(cos( 2.0 * M_PI * (x[0] / epsilon) ), 2.0) ) ) + stochastic perturbation
  // coeff.second = ( 0.1 + 1e-3 + ( 0.1 * sin( 2.0 * M_PI * (x[1] / epsilon) ) ) ) + stochastic perturbation
  const auto coeff = constants().coefficients_variant_A(x);

  flux[0][0] = coeff.first * direction_gradient[0][0] * (1.0 + pow(position_gradient[0][0], 2.0));
  flux[0][1] = coeff.second * direction_gradient[0][1] * (1.0 + pow(position_gradient[0][1], 2.0));
}

void ExactSolution::evaluate(const ExactSolution::DomainType&, ExactSolution::RangeType&) const {
  DUNE_THROW(Dune::NotImplemented, "Exact solution not available!");
}

void ExactSolution::jacobian(const ExactSolution::DomainType&, ExactSolution::JacobianRangeType&) const {
  DUNE_THROW(Dune::NotImplemented, "Exact solution not available!");
}

void ExactSolution::evaluate(const ExactSolution::DomainType& x, const TimeType&, ExactSolution::RangeType& y) const {
  evaluate(x, y);
}

} // namespace Five
} // namespace Problem
} // namespace Multiscale {
} // namespace Dune {
