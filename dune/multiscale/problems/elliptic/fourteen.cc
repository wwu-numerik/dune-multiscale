#include "fourteen.hh"

namespace Dune {
namespace Multiscale {
namespace Problem {
namespace Fourteen {

// default value for epsilon (if not sprecified in the parameter file)
CONSTANTSFUNCTION(0.05)

ModelProblemData::ModelProblemData() : IModelProblemData(constants()) {
  if (!constants().get("linear", true))
    DUNE_THROW(Dune::InvalidStateException, "problem Fourteen is entirely linear, but problem.linear was false");
  if (constants().get("stochastic_pertubation", false) && !(this->problemAllowsStochastics()))
    DUNE_THROW(Dune::InvalidStateException,
               "The problem does not allow stochastic perturbations. Please, switch the key off.");
}

std::string ModelProblemData::getMacroGridFile() const {
  return ("../dune/multiscale/grids/macro_grids/elliptic/cube_three.dgf");
}

bool ModelProblemData::problemIsPeriodic() const {
  return false; // = problem is periodic
}

bool ModelProblemData::problemAllowsStochastics() const {
  return false; // = problem does not allow stochastic perturbations
                // (if you want it, you must add the 'perturb' method provided
                // by 'constants.hh' - see model problems 4 to 7 for examples )
}

// evaluate f, i.e. return y=f(x) for a given x
// the following method defines 'f':
void FirstSource::evaluate(const DomainType& /*x*/, RangeType& y) const { y = 1.0; } // evaluate

void FirstSource::evaluate(const DomainType& x, const TimeType& /*time*/, RangeType& y) const { evaluate(x, y); }

Diffusion::Diffusion() {}

void Diffusion::diffusiveFlux(const DomainType& x, const JacobianRangeType& direction, JacobianRangeType& flux) const {

  const double eps = constants().epsilon;
  const double x0_eps = (x[0] / eps);
  const double cos_2_pi_x0_eps = cos(2.0 * M_PI * x0_eps);

  double coefficient = 1.1 + (0.5 * sin(std::floor(x[0] / eps)));
  coefficient += 0.5 * cos_2_pi_x0_eps;

  double a_0_0 = coefficient;
  double a_0_1 = 0.0;
  double a_1_0 = 0.0;
  double a_1_1 = coefficient;

  flux[0][0] = a_0_0 * direction[0][0] + a_0_1 * direction[0][1];
  flux[0][1] = a_1_0 * direction[0][0] + a_1_1 * direction[0][1];
} // diffusiveFlux

void Diffusion::jacobianDiffusiveFlux(const DomainType& x, const JacobianRangeType& /*position_gradient*/,
                                      const JacobianRangeType& direction_gradient, JacobianRangeType& flux) const {
  const double eps = constants().epsilon;
  const double x0_eps = (x[0] / eps);
  const double cos_2_pi_x0_eps = cos(2.0 * M_PI * x0_eps);

  double coefficient = 1.1 + (0.5 * sin(std::floor(x[0] / eps)));
  coefficient += 0.5 * cos_2_pi_x0_eps;

  double a_0_0 = coefficient;
  double a_0_1 = 0.0;
  double a_1_0 = 0.0;
  double a_1_1 = coefficient;

  flux[0][0] = a_0_0 * direction_gradient[0][0] + a_0_1 * direction_gradient[0][1];
  flux[0][1] = a_1_0 * direction_gradient[0][0] + a_1_1 * direction_gradient[0][1];
} // jacobianDiffusiveFlux

// evaluate Dirichlet Boundary Function
void DirichletBoundaryCondition::evaluate(const DomainType& x, RangeType& y) const {
  const double eps = constants().epsilon;
  const double x0_eps = (x[0] / eps);
  const double x1_eps = (x[1] / eps);
  const double sin_2_pi_x0_eps = sin(2.0 * M_PI * x0_eps);
  const double cos_2_pi_x1_eps = cos(2.0 * M_PI * x1_eps);

  y = sin_2_pi_x0_eps + cos_2_pi_x1_eps;
  y += 0.5 * exp(x[0] + x[1]);

} // evaluate

void DirichletBoundaryCondition::evaluate(const DomainType& x, const TimeType& /*time*/, RangeType& y) const {
  evaluate(x, y);
}

// evaluate Neumann Boundary Function 'q'
// q = A( \nabla u ) \cdot n    (on the Neumann boundary)
// ( A denotes the above diffusion operator and can be nonlinear )
void NeumannBoundaryCondition::evaluate(const DomainType& /*x*/, RangeType& y) const { y = 0.0; } // evaluate

void NeumannBoundaryCondition::evaluate(const DomainType& x, const TimeType& /*time*/, RangeType& y) const {
  evaluate(x, y);
}

void Dune::Multiscale::Problem::Fourteen::ExactSolution::evaluate(
    const Dune::Multiscale::Problem::Fourteen::ExactSolution::DomainType&,
    Dune::Multiscale::Problem::Fourteen::ExactSolution::RangeType&) const {
  DUNE_THROW(Dune::NotImplemented, "Exact solution not available!");
}

void Dune::Multiscale::Problem::Fourteen::ExactSolution::jacobian(
    const Dune::Multiscale::Problem::Fourteen::ExactSolution::DomainType&,
    Dune::Multiscale::Problem::Fourteen::ExactSolution::JacobianRangeType&) const {
  DUNE_THROW(Dune::NotImplemented, "Exact solution not available!");
}

void Dune::Multiscale::Problem::Fourteen::ExactSolution::evaluate(
    const Dune::Multiscale::Problem::Fourteen::ExactSolution::DomainType& x,
    const Dune::Multiscale::Problem::Fourteen::ExactSolution::TimeType&,
    Dune::Multiscale::Problem::Fourteen::ExactSolution::RangeType& y) const {
  evaluate(x, y);
}

} // namespace Fourteen
} // namespace Problem
} // namespace Multiscale {
} // namespace Dune {
