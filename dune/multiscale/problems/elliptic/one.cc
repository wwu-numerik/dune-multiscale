#include <config.h>
#include "one.hh"

namespace Dune {
namespace Multiscale {
namespace Problem {
namespace One {

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
  return false; // = problem does not allow stochastic perturbations
                // (if you want it, you must add the 'perturb' method provided
                // by 'constants.hh' - see model problems 4 to 7 for examples )
}

} // namespace One
} // namespace Problem
} // namespace Multiscale {
} // namespace Dune {

void Dune::Multiscale::Problem::One::Diffusion::diffusiveFlux(
    const Dune::Multiscale::Problem::One::Diffusion::DomainType& x,
    const Dune::Multiscale::Problem::One::Diffusion::JacobianRangeType& gradient,
    Dune::Multiscale::Problem::One::Diffusion::JacobianRangeType& flux) const {

  const auto diffusion_coefficient = 2.0 + sin(2.0 * M_PI * (x[0] / constants().epsilon));

  flux[0][0] = diffusion_coefficient * gradient[0][0];
  flux[0][1] = diffusion_coefficient * gradient[0][1];
}

void Dune::Multiscale::Problem::One::Diffusion::jacobianDiffusiveFlux(
    const Dune::Multiscale::Problem::One::Diffusion::DomainType& x,
    const Dune::Multiscale::Problem::One::Diffusion::JacobianRangeType&,
    const Dune::Multiscale::Problem::One::Diffusion::JacobianRangeType& direction_gradient,
    Dune::Multiscale::Problem::One::Diffusion::JacobianRangeType& flux) const {

  const auto diffusion_coefficient = 2.0 + sin(2.0 * M_PI * (x[0] / constants().epsilon));

  flux[0][0] = diffusion_coefficient * direction_gradient[0][0];
  flux[0][1] = diffusion_coefficient * direction_gradient[0][1];
}

void Dune::Multiscale::Problem::One::ExactSolution::evaluate(
    const Dune::Multiscale::Problem::One::ExactSolution::DomainType&,
    Dune::Multiscale::Problem::One::ExactSolution::RangeType&) const {
  DUNE_THROW(Dune::NotImplemented, "Exact solution not available!");
}

void Dune::Multiscale::Problem::One::ExactSolution::jacobian(
    const Dune::Multiscale::Problem::One::ExactSolution::DomainType&,
    Dune::Multiscale::Problem::One::ExactSolution::JacobianRangeType&) const {
  DUNE_THROW(Dune::NotImplemented, "Exact solution not available!");
}

void Dune::Multiscale::Problem::One::ExactSolution::evaluate(
    const Dune::Multiscale::Problem::One::ExactSolution::DomainType& x,
    const Dune::Multiscale::Problem::One::ExactSolution::TimeType&,
    Dune::Multiscale::Problem::One::ExactSolution::RangeType& y) const {
  evaluate(x, y);
}
