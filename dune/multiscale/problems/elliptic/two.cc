#include <config.h>
#include <assert.h>
#include <dune/common/exceptions.hh>
#include <dune/stuff/common/parameter/validation.hh>
#include <math.h>
#include <sstream>

#include "dune/multiscale/problems/base.hh"
#include "dune/multiscale/problems/constants.hh"
#include "two.hh"

namespace Dune {
namespace Multiscale {
namespace Problem {
namespace Two {

// default value for epsilon (if not sprecified in the parameter file)
CONSTANTSFUNCTION(0.05)

ModelProblemData::ModelProblemData() : IModelProblemData(constants()) {
  assert(constants_.epsilon != 0.0);
  if (constants().get("stochastic_pertubation", false) && !(this->problemAllowsStochastics()))
    DUNE_THROW(Dune::InvalidStateException,
               "The problem does not allow stochastic perturbations. Please, switch the key off.");
}

std::string ModelProblemData::getMacroGridFile() const {
  return ("../dune/multiscale/grids/macro_grids/elliptic/earth.dgf");
}

bool ModelProblemData::problemIsPeriodic() const {
  return false; // = problem is not periodic
}

bool ModelProblemData::problemAllowsStochastics() const { return false; }

void Diffusion::diffusiveFlux(
    const DomainType& x,
    const JacobianRangeType& gradient,
    JacobianRangeType& flux) const {

  double coefficient = 1.0 + (9.0 / 10.0) * sin(2.0 * M_PI * sqrt(fabs(2.0 * x[0])) / constants().epsilon) *
                                 sin(2.0 * M_PI * pow(1.5 * x[1], 2.0) / constants().epsilon);

  if (x[1] <= 0.3) {
    coefficient *= 4.0;
  }

  if ((x[1] > 0.3) && (x[1] < 0.6)) {
    coefficient *= 2.0 * (((-5.0 / 3.0) * x[1]) + (3.0 / 2.0));
  }

  if (x[1] >= 0.6) {
    coefficient *= 1.0;
  }

  flux[0][0] = coefficient * gradient[0][0];
  flux[0][1] = coefficient * gradient[0][1];
}

void Diffusion::jacobianDiffusiveFlux(
    const DomainType& x,
    const JacobianRangeType&,
    const JacobianRangeType& direction_gradient,
    JacobianRangeType& flux) const {

  double coefficient = 1.0 + (9.0 / 10.0) * sin(2.0 * M_PI * sqrt(fabs(2.0 * x[0])) / constants().epsilon) *
                                 sin(2.0 * M_PI * pow(1.5 * x[1], 2.0) / constants().epsilon);

  if (x[1] <= 0.3) {
    coefficient *= 4.0;
  }

  if ((x[1] > 0.3) && (x[1] < 0.6)) {
    coefficient *= 2.0 * (((-5.0 / 3.0) * x[1]) + (3.0 / 2.0));
  }

  if (x[1] >= 0.6) {
    coefficient *= 1.0;
  }

  flux[0][0] = coefficient * direction_gradient[0][0];
  flux[0][1] = coefficient * direction_gradient[0][1];
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

} // namespace two
} // namespace Problem
} // namespace Multiscale {
} // namespace Dune {
