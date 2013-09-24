#include "toy.hh"

namespace Dune {
namespace Multiscale {
namespace Problem {
namespace Toy {

// default value for epsilon (not required for this toy problem)
CONSTANTSFUNCTION(1.0)

ModelProblemData::ModelProblemData() : IModelProblemData(constants()) {
  if (!constants().get("linear", true))
    DUNE_THROW(Dune::InvalidStateException, "toy problem is entirely linear, but problem.linear was false");
  if (constants().get("stochastic_pertubation", false) && !(this->problemAllowsStochastics()))
    DUNE_THROW(Dune::InvalidStateException,
               "The problem does not allow stochastic perturbations. Please, switch the key off.");
}

std::string ModelProblemData::getMacroGridFile() const {
  return ("../dune/multiscale/grids/macro_grids/elliptic/cube_three.dgf");
}

bool ModelProblemData::problemIsPeriodic() const {
  return true; // = problem is periodic
}

bool ModelProblemData::problemAllowsStochastics() const {
  return false; // = problem does not allow stochastic perturbations
                // (if you want it, you must add the 'perturb' method provided
                // by 'constants.hh' - see model problems 4 to 7 for examples )
}

} // namespace Toy
} // namespace Problem
} // namespace Multiscale {
} // namespace Dune {

void Dune::Multiscale::Problem::Toy::ExactSolution::evaluate(
    const Dune::Multiscale::Problem::Toy::ExactSolution::DomainType& x,
    Dune::Multiscale::Problem::Toy::ExactSolution::RangeType& y) const {
  y = x[0] * (1.0 - x[0]) * (1.0 - x[1]) * x[1];
}

void Dune::Multiscale::Problem::Toy::ExactSolution::evaluate(
    const Dune::Multiscale::Problem::Toy::ExactSolution::DomainType& x,
    const Dune::Multiscale::Problem::Toy::ExactSolution::TimeType&,
    Dune::Multiscale::Problem::Toy::ExactSolution::RangeType& y) const {
  evaluate(x, y);
}

void Dune::Multiscale::Problem::Toy::ExactSolution::jacobian(
    const Dune::Multiscale::Problem::Toy::ExactSolution::DomainType& x,
    Dune::Multiscale::Problem::Toy::ExactSolution::JacobianRangeType& grad_u) const {
  grad_u[0][0] = (1.0 - x[0]) * (1.0 - x[1]) * x[1] - x[0] * (1.0 - x[1]) * x[1];
  grad_u[0][1] = x[0] * (1.0 - x[0]) * (1.0 - x[1]) - x[0] * (1.0 - x[0]) * x[1];
}

void Dune::Multiscale::Problem::Toy::Diffusion::diffusiveFlux(
    const Dune::Multiscale::Problem::Toy::Diffusion::DomainType& x,
    const Dune::Multiscale::Problem::Toy::Diffusion::JacobianRangeType& gradient,
    Dune::Multiscale::Problem::Toy::Diffusion::JacobianRangeType& flux) const {
  double a_0 = 1.0 + pow(x[0], 2.0);

  flux[0][0] = a_0 * gradient[0][0];
  flux[0][1] = a_0 * gradient[0][1];
}
