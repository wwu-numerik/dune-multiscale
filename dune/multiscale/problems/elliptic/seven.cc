#include "seven.hh"

namespace Dune {
namespace Multiscale {
namespace Problem {
namespace Seven {

// default value for epsilon (if not sprecified in the parameter file)
CONSTANTSFUNCTION(0.05)

ModelProblemData::ModelProblemData() : IModelProblemData(constants()) {
  if (constants().get("linear", true))
    DUNE_THROW(Dune::InvalidStateException, "Problem seven is entirely nonlinear, but problem.linear was true.");
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

} // namespace Seven
} // namespace Problem
} // namespace Multiscale {
} // namespace Dune {

void Dune::Multiscale::Problem::Seven::Diffusion::diffusiveFlux(
    const Dune::Multiscale::Problem::Seven::Diffusion::DomainType& x,
    const Dune::Multiscale::Problem::Seven::Diffusion::JacobianRangeType& gradient,
    Dune::Multiscale::Problem::Seven::Diffusion::JacobianRangeType& flux) const {
  const auto coeff = constants().coefficients(x);
  flux[0][0] = coeff.first * (gradient[0][0] + ((1.0 / 3.0) * pow(gradient[0][0], 3.0)));
  flux[0][1] = coeff.second * (gradient[0][1] + ((1.0 / 3.0) * pow(gradient[0][1], 3.0)));
}

void Dune::Multiscale::Problem::Seven::Diffusion::jacobianDiffusiveFlux(
    const Dune::Multiscale::Problem::Seven::Diffusion::DomainType& x,
    const Dune::Multiscale::Problem::Seven::Diffusion::JacobianRangeType& position_gradient,
    const Dune::Multiscale::Problem::Seven::Diffusion::JacobianRangeType& direction_gradient,
    Dune::Multiscale::Problem::Seven::Diffusion::JacobianRangeType& flux) const {
  const auto coeff = constants().coefficients(x);
  flux[0][0] = coeff.first * direction_gradient[0][0] * (1.0 + pow(position_gradient[0][0], 2.0));
  flux[0][1] = coeff.second * direction_gradient[0][1] * (1.0 + pow(position_gradient[0][1], 2.0));
}

void Dune::Multiscale::Problem::Seven::ExactSolution::evaluate(
    const Dune::Multiscale::Problem::Seven::ExactSolution::DomainType&,
    Dune::Multiscale::Problem::Seven::ExactSolution::RangeType&) const {
  DUNE_THROW(Dune::NotImplemented, "Exact solution not available!");
}

void Dune::Multiscale::Problem::Seven::ExactSolution::jacobian(
    const Dune::Multiscale::Problem::Seven::ExactSolution::DomainType&,
    Dune::Multiscale::Problem::Seven::ExactSolution::JacobianRangeType&) const {
  DUNE_THROW(Dune::NotImplemented, "Exact solution not available!");
}

void Dune::Multiscale::Problem::Seven::ExactSolution::evaluate(
    const Dune::Multiscale::Problem::Seven::ExactSolution::DomainType& x,
    const Dune::Multiscale::Problem::Seven::ExactSolution::TimeType&,
    Dune::Multiscale::Problem::Seven::ExactSolution::RangeType& y) const {
  evaluate(x, y);
}
