#include <config.h>
#include <dune/common/exceptions.hh>
#include <dune/stuff/common/parameter/validation.hh>
#include <sstream>

#include "dune/multiscale/problems/base.hh"
#include "dune/multiscale/problems/constants.hh"
#include "toy.hh"

namespace Dune {
namespace Multiscale {
namespace Problem {
namespace Toy {

// default value for epsilon (not required for this toy problem)
CONSTANTSFUNCTION(1.0)

ModelProblemData::ModelProblemData() : IModelProblemData(constants()) {
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

bool ModelProblemData::problemAllowsStochastics() const { return false; }

void FirstSource::evaluate(const DomainType& x, RangeType& y) const {
  double a_0_x_0 = 1.0 + pow(x[0], 2.0);
  double a_1_x_1 = 1.0 + pow(x[0], 2.0);

  double grad_a_0_x_0 = 2.0 * x[0];
  double grad_a_1_x_1 = 0.0;

  typename FunctionSpaceType::JacobianRangeType grad_u(0.0);

  grad_u[0][0] = (1.0 - x[0]) * (1.0 - x[1]) * x[1] - x[0] * (1.0 - x[1]) * x[1];
  grad_u[0][1] = x[0] * (1.0 - x[0]) * (1.0 - x[1]) - x[0] * (1.0 - x[0]) * x[1];

  const typename FunctionSpaceType::RangeType d_xx_u = (-2.0) * (1.0 - x[1]) * x[1];
  const typename FunctionSpaceType::RangeType d_yy_u = (-2.0) * (1.0 - x[0]) * x[0];

  y = 0.0;
  y -= grad_a_0_x_0 * grad_u[0][0];
  y -= a_0_x_0 * d_xx_u;

  y -= grad_a_1_x_1 * grad_u[0][1];
  y -= a_1_x_1 * d_yy_u;
}

void ExactSolution::evaluate(const ExactSolution::DomainType& x, ExactSolution::RangeType& y) const {
  y = x[0] * (1.0 - x[0]) * (1.0 - x[1]) * x[1];
}

void ExactSolution::evaluate(const ExactSolution::DomainType& x, const TimeType&, ExactSolution::RangeType& y) const {
  evaluate(x, y);
}

void ExactSolution::jacobian(const ExactSolution::DomainType& x, ExactSolution::JacobianRangeType& grad_u) const {
  grad_u[0][0] = (1.0 - x[0]) * (1.0 - x[1]) * x[1] - x[0] * (1.0 - x[1]) * x[1];
  grad_u[0][1] = x[0] * (1.0 - x[0]) * (1.0 - x[1]) - x[0] * (1.0 - x[0]) * x[1];
}

void Diffusion::diffusiveFlux(const DomainType& x, const JacobianRangeType& gradient, JacobianRangeType& flux) const {
  double a_0 = 1.0 + pow(x[0], 2.0);

  flux[0][0] = a_0 * gradient[0][0];
  flux[0][1] = a_0 * gradient[0][1];
}

void DirichletData::evaluate(const DomainType& x, RangeType& y) const {
  y = x[0] * (1.0 - x[0]) * (1.0 - x[1]) * x[1];
}

void DirichletData::evaluate(const DomainType& x, const TimeType&, RangeType& y) const {
  evaluate(x, y);
}

void DirichletData::jacobian(const DomainType& x, JacobianRangeType& grad_u) const {
  grad_u[0][0] = (1.0 - x[0]) * (1.0 - x[1]) * x[1] - x[0] * (1.0 - x[1]) * x[1];
  grad_u[0][1] = x[0] * (1.0 - x[0]) * (1.0 - x[1]) - x[0] * (1.0 - x[0]) * x[1];
}

} // namespace Toy
} // namespace Problem
} // namespace Multiscale {
} // namespace Dune {
