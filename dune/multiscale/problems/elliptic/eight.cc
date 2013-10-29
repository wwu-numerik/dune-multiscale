#include <config.h>
#include <dune/common/exceptions.hh>
#include <dune/stuff/common/parameter/validation.hh>
#include <math.h>
#include <sstream>

#include "dune/multiscale/problems/base.hh"
#include "dune/multiscale/problems/constants.hh"
#include "eight.hh"

namespace Dune {
namespace Multiscale {
namespace Problem {
namespace Eight {

// default value for epsilon
CONSTANTSFUNCTION(0.001)

ModelProblemData::ModelProblemData() : IModelProblemData(constants()) {
  if (constants().get("linear", true))
    DUNE_THROW(Dune::InvalidStateException, "problem eight is entirely nonlinear, but problem.linear was true");
  if (constants().get("stochastic_pertubation", false) && !(this->problemAllowsStochastics()))
    DUNE_THROW(Dune::InvalidStateException,
               "The problem does not allow stochastic perturbations. Please, switch the key off.");
}

std::string ModelProblemData::getMacroGridFile() const {
  return ("../dune/multiscale/grids/macro_grids/elliptic/unit_cube.dgf");
}

bool ModelProblemData::problemIsPeriodic() const {
  return false; // = problem is not periodic
}

bool ModelProblemData::problemAllowsStochastics() const {
  return false; // = problem does not allow stochastic perturbations
                // (if you want it, you must add the 'perturb' method provided
                // by 'constants.hh' - see model problems 4 to 7 for examples )
}

void ExactSolution::evaluate(const ExactSolution::DomainType& x, ExactSolution::RangeType& y) const {
  y = (-1.0) * ((x[0] * x[0]) - x[0]) * ((x[1] * x[1]) - x[1]);
  y -= constants().epsilon * (x[0] + x[1]) * sin(2.0 * M_PI * x[0] / constants().epsilon) *
       sin(2.0 * M_PI * x[1] / constants().epsilon);
}

void ExactSolution::evaluate(const ExactSolution::DomainType& x, const ExactSolution::TimeType&,
                             ExactSolution::RangeType& y) const {
  evaluate(x, y);
}

void ExactSolution::jacobian(const ExactSolution::DomainType&,
                             ExactSolution::FunctionSpaceType::JacobianRangeType&) const {
  DUNE_THROW(Dune::NotImplemented, "Dummy body for all-problem compile");
}

void Diffusion::jacobianDiffusiveFlux(const Diffusion::DomainType& x,
                                      const Diffusion::JacobianRangeType& position_gradient,
                                      const Diffusion::JacobianRangeType& direction_gradient,
                                      Diffusion::JacobianRangeType& flux) const {
  double coefficient = 2.0 + sin(2.0 * M_PI * (x[0] + x[1]) / constants().epsilon);

  flux[0][0] = direction_gradient[0][0] * (1.0 + 3.0 * coefficient * pow(position_gradient[0][0], 2.0));
  flux[0][1] = direction_gradient[0][1] * (1.0 + 3.0 * coefficient * pow(position_gradient[0][1], 2.0));

  flux[0][0] *= (-1.0);
  flux[0][1] *= (-1.0);
}

void Diffusion::diffusiveFlux(const Diffusion::DomainType& x, const Diffusion::JacobianRangeType& gradient,
                              Diffusion::JacobianRangeType& flux) const {
  double coefficient = 2.0 + sin(2.0 * M_PI * (x[0] + x[1]) / constants().epsilon);

  flux[0][0] = gradient[0][0] + (coefficient * pow(gradient[0][0], 3.0));
  flux[0][0] -= additivePart(x, 0, 1);
  flux[0][0] *= (-1.0);

  flux[0][1] = gradient[0][1] + (coefficient * pow(gradient[0][1], 3.0));
  flux[0][1] -= additivePart(x, 1, 0);
  flux[0][1] *= (-1.0);
}

double Diffusion::additivePart(const Diffusion::DomainType& x, const int i, const int j) const {
  double y = 0.0;

  y -= (x[0] + x[1]) * cos(2.0 * M_PI * x[i] / constants().epsilon) * sin(2.0 * M_PI * x[j] / constants().epsilon);

  double helper1 = 1.0;
  helper1 *= (2.0 + sin(2.0 * M_PI * (x[0] + x[1]) / constants().epsilon));

  double helper2 = 1.0;
  helper2 *= 3.0;
  helper2 *=
      ((2.0 * x[i] - 1.0) * (x[j] * x[j] - x[j])) +
      ((x[0] + x[1]) * cos(2.0 * M_PI * x[i] / constants().epsilon) * sin(2.0 * M_PI * x[j] / constants().epsilon));
  helper2 *= (2.0 * x[i] - 1.0) * (x[j] * x[j] - x[j]) * (x[0] + x[1]) * cos(2.0 * M_PI * x[i] / constants().epsilon) *
             sin(2.0 * M_PI * x[j] / constants().epsilon);
  helper2 += pow(
      (x[0] + x[1]) * cos(2.0 * M_PI * x[i] / constants().epsilon) * sin(2.0 * M_PI * x[j] / constants().epsilon), 3.0);

  helper1 *= helper2;
  y -= helper1;
  y -= sin(2.0 * M_PI * (x[0] + x[1]) / constants().epsilon) * pow((2.0 * x[i]) - 1.0, 3.0) *
       pow((x[j] * x[j]) - x[j], 3.0);
  return y;
}

void FirstSource::evaluate(const FirstSource::DomainType& x, FirstSource::RangeType& y) const {
  y = 0.0;
  y += 2.0 * (x[0] + x[1] - pow(x[0], 2.0) - pow(x[1], 2.0));
  y -= 12.0 * pow((2.0 * x[0]) - 1.0, 2.0) * pow((x[1] * x[1]) - x[1], 3.0);
  y -= 12.0 * pow((2.0 * x[1]) - 1.0, 2.0) * pow((x[0] * x[0]) - x[0], 3.0);
}

void FirstSource::evaluate(const FirstSource::DomainType& x, const FirstSource::TimeType&,
                           FirstSource::RangeType& y) const {
  evaluate(x, y);
}

} // namespace Eight
} // namespace Problem
} // namespace Multiscale {
} // namespace Dune {
