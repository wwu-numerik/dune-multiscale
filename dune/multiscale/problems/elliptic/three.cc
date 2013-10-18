#include <config.h>
#include "three.hh"

namespace Dune {
namespace Multiscale {
namespace Problem {
namespace Three {

// default value for epsilon (if not sprecified in the parameter file)
CONSTANTSFUNCTION(0.05)

ModelProblemData::ModelProblemData() : IModelProblemData(constants()) {
  assert(constants_.epsilon != 0.0);
  if (constants().get("linear", true))
    DUNE_THROW(Dune::InvalidStateException, "problem three is entirely nonlinear, but problem.linear was true");
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

} // namespace Three
} // namespace Problem
} // namespace Multiscale {
} // namespace Dune {

void Dune::Multiscale::Problem::Three::Diffusion::diffusiveFlux(
    const Dune::Multiscale::Problem::Three::Diffusion::DomainType& x,
    const Dune::Multiscale::Problem::Three::Diffusion::JacobianRangeType& gradient,
    Dune::Multiscale::Problem::Three::Diffusion::JacobianRangeType& flux) const {
  double coefficient = 1.0 + (9.0 / 10.0) * sin(2.0 * M_PI * sqrt(fabs(2.0 * x[0])) / constants().epsilon) *
                                 sin(2.0 * M_PI * pow(1.5 * x[1], 2.0) / constants().epsilon);

  if ((x[1] > 0.3) && (x[1] < 0.6))
    coefficient *= ((-3.0) * x[1] + 1.9);

  if (x[1] >= 0.6)
    coefficient *= 0.1;

  flux[0][0] = coefficient * (gradient[0][0] + ((1.0 / 3.0) * pow(gradient[0][0], 3.0)));
  flux[0][1] = coefficient * (gradient[0][1] + ((1.0 / 3.0) * pow(gradient[0][1], 3.0)));
}

void Dune::Multiscale::Problem::Three::Diffusion::jacobianDiffusiveFlux(
    const Dune::Multiscale::Problem::Three::Diffusion::DomainType& x,
    const Dune::Multiscale::Problem::Three::Diffusion::JacobianRangeType& position_gradient,
    const Dune::Multiscale::Problem::Three::Diffusion::JacobianRangeType& direction_gradient,
    Dune::Multiscale::Problem::Three::Diffusion::JacobianRangeType& flux) const {
  double coefficient = 1.0 + (9.0 / 10.0) * sin(2.0 * M_PI * sqrt(fabs(2.0 * x[0])) / constants().epsilon) *
                                 sin(2.0 * M_PI * pow(1.5 * x[1], 2.0) / constants().epsilon);

  if ((x[1] > 0.3) && (x[1] < 0.6)) {
    coefficient *= ((-3.0) * x[1] + 1.9);
  }

  if (x[1] >= 0.6) {
    coefficient *= 0.1;
  }

  flux[0][0] = coefficient * direction_gradient[0][0] * (1.0 + pow(position_gradient[0][0], 2.0));
  flux[0][1] = coefficient * direction_gradient[0][1] * (1.0 + pow(position_gradient[0][1], 2.0));
}

void Dune::Multiscale::Problem::Three::ExactSolution::evaluate(
    const Dune::Multiscale::Problem::Three::ExactSolution::DomainType&,
    Dune::Multiscale::Problem::Three::ExactSolution::RangeType&) const {
  DUNE_THROW(Dune::NotImplemented, "Exact solution not available!");
}

void Dune::Multiscale::Problem::Three::ExactSolution::jacobian(
    const Dune::Multiscale::Problem::Three::ExactSolution::DomainType&,
    Dune::Multiscale::Problem::Three::ExactSolution::JacobianRangeType&) const {
  DUNE_THROW(Dune::NotImplemented, "Exact solution not available!");
}

void Dune::Multiscale::Problem::Three::ExactSolution::evaluate(
    const Dune::Multiscale::Problem::Three::ExactSolution::DomainType& x,
    const Dune::Multiscale::Problem::Three::ExactSolution::TimeType&,
    Dune::Multiscale::Problem::Three::ExactSolution::RangeType& y) const {
  evaluate(x, y);
}
