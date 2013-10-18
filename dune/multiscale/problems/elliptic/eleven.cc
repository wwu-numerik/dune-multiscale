#include <config.h>
#include "eleven.hh"

namespace Dune {
namespace Multiscale {
namespace Problem {
namespace Eleven {

// default value for epsilon (if not sprecified in the parameter file)
CONSTANTSFUNCTION(0.05)

ModelProblemData::ModelProblemData() : IModelProblemData(constants()) {
  if (constants().get("linear", true))
    DUNE_THROW(Dune::InvalidStateException, "problem eleven is entirely nonlinear, but problem.linear was true");
  if (constants().get("stochastic_pertubation", false) && !(this->problemAllowsStochastics()))
    DUNE_THROW(Dune::InvalidStateException,
               "The problem does not allow stochastic perturbations. Please, switch the key off.");
}

std::string ModelProblemData::getMacroGridFile() const {
  return ("../dune/multiscale/grids/macro_grids/elliptic/cube_three_dirichlet_neumann.dgf");
}

bool ModelProblemData::problemIsPeriodic() const {
  return true; // = problem is periodic
}

bool ModelProblemData::problemAllowsStochastics() const {
  return false; // = problem does not allow stochastic perturbations
                // (if you want it, you must add the 'perturb' method provided
                // by 'constants.hh' - see model problems 4 to 7 for examples )
}

// evaluate f, i.e. return y=f(x) for a given x
// the following method defines 'f':
void FirstSource::evaluate(const DomainType& x, RangeType& y) const {
  y = -6.0 * exp(x[0] + x[1]);

  if (exp(x[0] + x[1]) <= (-1.0)) {
    y += (-1.0) * (1.0 + pow(x[0] + x[1], 2.0));
  } else if (exp(x[0] + x[1]) >= 1) {
    y += (1.0 + pow(x[0] + x[1], 2.0));
  } else {
    y += (1.0 + pow(x[0] + x[1], 2.0)) * sin(0.5 * M_PI * exp(x[0] + x[1]));
  }

} // evaluate

void FirstSource::evaluate(const DomainType& x, const TimeType& /*time*/, RangeType& y) const { evaluate(x, y); }

Diffusion::Diffusion() {}

void Diffusion::diffusiveFlux(const DomainType& /*x*/, const JacobianRangeType& direction,
                              JacobianRangeType& flux) const {
  double a_0_0 = 3.0;
  double a_0_1 = 0.0;
  double a_1_0 = 0.0;
  double a_1_1 = 3.0;

  flux[0][0] = a_0_0 * direction[0][0] + a_0_1 * direction[0][1];
  flux[0][1] = a_1_0 * direction[0][0] + a_1_1 * direction[0][1];
} // diffusiveFlux

void Diffusion::jacobianDiffusiveFlux(const DomainType& /*x*/, const JacobianRangeType& /*position_gradient*/,
                                      const JacobianRangeType& direction_gradient, JacobianRangeType& flux) const {
  double a_0_0 = 3.0;
  double a_0_1 = 0.0;
  double a_1_0 = 0.0;
  double a_1_1 = 3.0;

  flux[0][0] = a_0_0 * direction_gradient[0][0] + a_0_1 * direction_gradient[0][1];
  flux[0][1] = a_1_0 * direction_gradient[0][0] + a_1_1 * direction_gradient[0][1];
} // jacobianDiffusiveFlux

// dummy
void LowerOrderTerm::evaluate(const DomainType& /*x*/, RangeType& /*y*/) const {
  DUNE_THROW(Dune::NotImplemented, "Inadmissible call for 'evaluate'");
}

// dummy
void LowerOrderTerm::evaluate(const DomainType& /*x*/, const TimeType& /*time*/, RangeType& /*y*/) const {
  DUNE_THROW(Dune::NotImplemented, "Inadmissible call for 'evaluate'");
}

// evaluate y = F(x, u(x), \grad u(x))
// 'position = u(x)', 'direction_gradient = \grad u(x)'
//! rename this into 'lower order term' or something similar
void LowerOrderTerm::evaluate(const DomainType& x, const RangeType& position,
                              const JacobianRangeType& /*direction_gradient*/, RangeType& y) const {

  // F(x,p, (z_1,z_2) ) = c g(p) z_2
  y = 0.0;

  if (position <= (-1.0)) {
    y = (-1.0) * (1.0 + pow(x[0] + x[1], 2.0));
  } else if (position >= 1) {
    y = (1.0 + pow(x[0] + x[1], 2.0));
  } else {
    y = (1.0 + pow(x[0] + x[1], 2.0)) * sin(0.5 * M_PI * position);
  }

} // evaluate

// evaluate position derivative y = d_1 F (x, u(x), \grad u(x))  (derivative with respect to the second componenent
// 'u(x)')
// 'position = u(x)', 'direction_gradient = \grad u(x)'
void LowerOrderTerm::position_derivative(const DomainType& x, const RangeType& position,
                                         const JacobianRangeType& /*direction_gradient*/, RangeType& y) const {

  // \partial_p F(x,p, (z_1,z_2) ) = ..
  y = 0.0;

  if (position <= (-1.0)) {
    y = 0.0;
  } else if (position >= 1) {
    y = 0.0;
  } else {
    y = 0.5 * M_PI * (1.0 + pow(x[0] + x[1], 2.0)) * cos(0.5 * M_PI * position);
  }

} // position_derivative

// evaluate position derivative y = d_2 F (x, u(x), \grad u(x))  (derivative with respect to the third componenent 'grad
// u(x)')
// 'position = u(x)', 'direction_gradient = \grad u(x)'
void LowerOrderTerm::direction_derivative(const DomainType& /*x*/, const RangeType& /*position*/,
                                          const JacobianRangeType& /*direction_gradient*/, JacobianRangeType& y) const {
  // \grad_z F(x,p, (z_1,z_2) ) = ...
  y[0][1] = 0.0;
  y[0][0] = 0.0;

} // direction_derivative

// evaluate Dirichlet Boundary Function
void DirichletBoundaryCondition::evaluate(const DomainType& x, RangeType& y) const { y = exp(x[0] + x[1]); } // evaluate

void DirichletBoundaryCondition::evaluate(const DomainType& x, const TimeType& /*time*/, RangeType& y) const {
  evaluate(x, y);
}

// evaluate Neumann Boundary Function 'q'
// q = A( \nabla u ) \cdot n    (on the Neumann boundary)
// ( A denotes the above diffusion operator and can be nonlinear )
void NeumannBoundaryCondition::evaluate(const DomainType& x, RangeType& y) const {
  y = (-3.0) * exp(x[0] + x[1]);

} // evaluate

void NeumannBoundaryCondition::evaluate(const DomainType& x, const TimeType& /*time*/, RangeType& y) const {
  evaluate(x, y);
}

ExactSolution::ExactSolution() {}

void ExactSolution::evaluate(const DomainType& x, RangeType& y) const { y = exp(x[0] + x[1]); } // evaluate

void ExactSolution::jacobian(const DomainType& x, typename FunctionSpaceType::JacobianRangeType& grad_u) const {

  grad_u[0][0] = exp(x[0] + x[1]);
  grad_u[0][1] = exp(x[0] + x[1]);

} // jacobian

void ExactSolution::evaluate(const DomainType& x, const TimeType& /*timedummy*/, RangeType& y) const { evaluate(x, y); }
} // namespace Eleven
} // namespace Problem
} // namespace Multiscale {
} // namespace Dune {
