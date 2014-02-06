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

bool ModelProblemData::problemAllowsStochastics() const {
  return false; // = problem does not allow stochastic perturbations
                // (if you want it, you must add the 'perturb' method provided
                // by 'constants.hh' - see model problems 4 to 7 for examples )
}

} // namespace two
} // namespace Problem
} // namespace Multiscale {
} // namespace Dune {

void Dune::Multiscale::Problem::Two::Diffusion::diffusiveFlux(
    const Dune::Multiscale::Problem::Two::Diffusion::DomainType& x,
    const Dune::Multiscale::Problem::Two::Diffusion::JacobianRangeType& gradient,
    Dune::Multiscale::Problem::Two::Diffusion::JacobianRangeType& flux) const {

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

void Dune::Multiscale::Problem::Two::Diffusion::jacobianDiffusiveFlux(
    const Dune::Multiscale::Problem::Two::Diffusion::DomainType& x,
    const Dune::Multiscale::Problem::Two::Diffusion::JacobianRangeType&,
    const Dune::Multiscale::Problem::Two::Diffusion::JacobianRangeType& direction_gradient,
    Dune::Multiscale::Problem::Two::Diffusion::JacobianRangeType& flux) const {

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

void Dune::Multiscale::Problem::Two::ExactSolution::evaluate(
    const Dune::Multiscale::Problem::Two::ExactSolution::DomainType&,
    Dune::Multiscale::Problem::Two::ExactSolution::RangeType&) const {
  DUNE_THROW(Dune::NotImplemented, "Exact solution not available!");
}

void Dune::Multiscale::Problem::Two::ExactSolution::jacobian(
    const Dune::Multiscale::Problem::Two::ExactSolution::DomainType&,
    Dune::Multiscale::Problem::Two::ExactSolution::JacobianRangeType&) const {
  DUNE_THROW(Dune::NotImplemented, "Exact solution not available!");
}

void Dune::Multiscale::Problem::Two::ExactSolution::evaluate(
    const Dune::Multiscale::Problem::Two::ExactSolution::DomainType& x,
    const TimeType&,
    Dune::Multiscale::Problem::Two::ExactSolution::RangeType& y) const {
  evaluate(x, y);
}
