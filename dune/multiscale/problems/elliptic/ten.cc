#include <config.h>
#include <assert.h>
#include <dune/common/exceptions.hh>
#include <dune/stuff/common/parameter/validation.hh>
#include <math.h>
#include <sstream>

#include "dune/multiscale/problems/base.hh"
#include "dune/multiscale/problems/constants.hh"
#include "ten.hh"

namespace Dune {
namespace Multiscale {
namespace Problem {
namespace Ten {

// description see below 0.05
CONSTANTSFUNCTION(0.05)

ModelProblemData::ModelProblemData() : IModelProblemData(constants()) {
  assert(constants_.epsilon == 0.0);
  if (constants().get("stochastic_pertubation", false) && !(this->problemAllowsStochastics()))
    DUNE_THROW(Dune::InvalidStateException,
               "The problem does not allow stochastic perturbations. Please, switch the key off.");
}

std::string ModelProblemData::getMacroGridFile() const {
  return ("../dune/multiscale/grids/macro_grids/elliptic/cube_three.dgf"); // _strange_grid
}

bool ModelProblemData::problemIsPeriodic() const {
  return false; // = problem is not periodic
}

bool ModelProblemData::problemAllowsStochastics() const {
  return false; // = problem does not allow stochastic perturbations
                // (if you want it, you must add the 'perturb' method provided
                // by 'constants.hh' - see model problems 4 to 7 for examples )
}

} // namespace Ten
} // namespace Problem
} // namespace Multiscale {
} // namespace Dune {

void Dune::Multiscale::Problem::Ten::Diffusion::diffusiveFlux(
    const Dune::Multiscale::Problem::Ten::Diffusion::DomainType& x,
    const Dune::Multiscale::Problem::Ten::Diffusion::JacobianRangeType& gradient,
    Dune::Multiscale::Problem::Ten::Diffusion::JacobianRangeType& flux) const {
  double diff_coef = 0.0;

  double coefficient =
      (1.0 / (8.0 * M_PI * M_PI)) *
      (1.0 + (0.5 * cos(2.0 * M_PI * (x[0] / constants().epsilon)) * sin(2.0 * M_PI * (x[1] / constants().epsilon))));

  double constant_val = 0.0005;

  double r1 = 0.425;
  double r2 = 0.125;

  // check if x part of the ellipse
  double position = (pow(x[0] - 0.5, 2.0) / pow(r1, 2.0)) + (pow(x[1] - 0.5, 2.0) / pow(r2, 2.0));

  if (position <= 1.0) {
    if (position >= 0.9) {
      diff_coef = ((10.0 * position) - 9.0) * coefficient;
      diff_coef += (1.0 - ((10.0 * position) - 9.0)) * constant_val;
    } else {
      diff_coef = constant_val;

      r1 = 0.35;
      r2 = 0.022;

      double new_position = (pow(x[0] - 0.5, 2.0) / pow(r1, 2.0)) + (pow(x[1] - 0.5, 2.0) / pow(r2, 2.0));

      if (new_position <= 1.0)
        diff_coef *= 100.0;

      r1 = 0.25;
      r2 = 0.022;

      new_position = (pow(x[0] - 0.5, 2.0) / pow(r1, 2.0)) + (pow(x[1] - 0.56, 2.0) / pow(r2, 2.0));

      if (new_position <= 1.0)
        diff_coef *= 100.0;

      new_position = (pow(x[0] - 0.5, 2.0) / pow(r1, 2.0)) + (pow(x[1] - 0.44, 2.0) / pow(r2, 2.0));

      if (new_position <= 1.0)
        diff_coef *= 100.0;
    }
  } else {
    diff_coef = coefficient;
  }

  const double stab = 0.0;
  flux[0][0] = diff_coef * gradient[0][0] + stab * gradient[0][1];
  flux[0][1] = diff_coef * gradient[0][1] + stab * gradient[0][0];
}

void Dune::Multiscale::Problem::Ten::Diffusion::jacobianDiffusiveFlux(
    const Dune::Multiscale::Problem::Ten::Diffusion::DomainType& x,
    const Dune::Multiscale::Problem::Ten::Diffusion::JacobianRangeType&,
    const Dune::Multiscale::Problem::Ten::Diffusion::JacobianRangeType& direction_gradient,
    Dune::Multiscale::Problem::Ten::Diffusion::JacobianRangeType& flux) const {
  double diff_coef = 0.0;

  double coefficient =
      (1.0 / (8.0 * M_PI * M_PI)) *
      (1.0 + (0.5 * cos(2.0 * M_PI * (x[0] / constants().epsilon)) * sin(2.0 * M_PI * (x[1] / constants().epsilon))));

  double constant_val = 0.0005;

  double r1 = 0.425;
  double r2 = 0.125;

  // check if x part of the ellipse
  double position = (pow(x[0] - 0.5, 2.0) / pow(r1, 2.0)) + (pow(x[1] - 0.5, 2.0) / pow(r2, 2.0));

  if (position <= 1.0) {
    if (position >= 0.9) {
      diff_coef = ((10.0 * position) - 9.0) * coefficient;
      diff_coef += (1.0 - ((10.0 * position) - 9.0)) * constant_val;
    } else {
      diff_coef = constant_val;

      r1 = 0.35;
      r2 = 0.022;

      double new_position = (pow(x[0] - 0.5, 2.0) / pow(r1, 2.0)) + (pow(x[1] - 0.5, 2.0) / pow(r2, 2.0));

      if (new_position <= 1.0)
        diff_coef *= 100.0;

      r1 = 0.25;
      r2 = 0.022;

      new_position = (pow(x[0] - 0.5, 2.0) / pow(r1, 2.0)) + (pow(x[1] - 0.56, 2.0) / pow(r2, 2.0));

      if (new_position <= 1.0)
        diff_coef *= 100.0;

      new_position = (pow(x[0] - 0.5, 2.0) / pow(r1, 2.0)) + (pow(x[1] - 0.44, 2.0) / pow(r2, 2.0));

      if (new_position <= 1.0)
        diff_coef *= 100.0;
    }
  } else {
    diff_coef = coefficient;
  }
  flux[0][0] = diff_coef * direction_gradient[0][0];
  flux[0][1] = diff_coef * direction_gradient[0][1];
}

void Dune::Multiscale::Problem::Ten::ExactSolution::evaluate(
    const Dune::Multiscale::Problem::Ten::ExactSolution::DomainType&,
    Dune::Multiscale::Problem::Ten::ExactSolution::RangeType&) const {
  DUNE_THROW(Dune::NotImplemented, "Exact solution not available!");
}

void Dune::Multiscale::Problem::Ten::ExactSolution::jacobian(
    const Dune::Multiscale::Problem::Ten::ExactSolution::DomainType&,
    Dune::Multiscale::Problem::Ten::ExactSolution::JacobianRangeType&) const {
  DUNE_THROW(Dune::NotImplemented, "Exact solution not available!");
}

void Dune::Multiscale::Problem::Ten::ExactSolution::evaluate(
    const Dune::Multiscale::Problem::Ten::ExactSolution::DomainType& x,
    const TimeType&,
    Dune::Multiscale::Problem::Ten::ExactSolution::RangeType& y) const {
  evaluate(x, y);
}
