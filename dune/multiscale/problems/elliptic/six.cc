#include "six.hh"

namespace Dune {
namespace Multiscale {
namespace Problem {
namespace Six {

// default value for epsilon (if not sprecified in the parameter file)
CONSTANTSFUNCTION( 0.05 )

ModelProblemData::ModelProblemData()
  : IModelProblemData(constants()) {
    assert( constants_.epsilon != 0.0);
    if (constants().get("stochastic_pertubation", false) && !(this->problemAllowsStochastics()) )
       DUNE_THROW(Dune::InvalidStateException, "The problem does not allow stochastic perturbations. Please, switch the key off.");
}

std::string ModelProblemData::getMacroGridFile() const {
  return("../dune/multiscale/grids/macro_grids/elliptic/cube_two.dgf");
}

bool ModelProblemData::problemIsPeriodic() const {
  return true; // = problem is periodic
}

bool ModelProblemData::problemAllowsStochastics() const {
  return true; // = problem allows stochastic perturbations
}


} //namespace Six
} //namespace Problem
} //namespace Multiscale {
} //namespace Dune {

void Dune::Multiscale::Problem::Six::Diffusion::diffusiveFlux(const Dune::Multiscale::Problem::Six::Diffusion::DomainType &x, const Dune::Multiscale::Problem::Six::Diffusion::JacobianRangeType &gradient, Dune::Multiscale::Problem::Six::Diffusion::JacobianRangeType &flux) const {

  // coeff.first = 1.01 + cos( 2.0 * M_PI * (x[0] / epsilon) ) + stochastic perturbation
  // coeff.second = 1.01 + cos( 2.0 * M_PI * (x[0] / epsilon) ) + stochastic perturbation
  const auto coeff = constants().coefficients(x);
  flux[0][0] = coeff.first * gradient[0][0];
  flux[0][1] = coeff.second * gradient[0][1];
}

void Dune::Multiscale::Problem::Six::Diffusion::jacobianDiffusiveFlux(const Dune::Multiscale::Problem::Six::Diffusion::DomainType &x, const Dune::Multiscale::Problem::Six::Diffusion::JacobianRangeType &, const Dune::Multiscale::Problem::Six::Diffusion::JacobianRangeType &direction_gradient, Dune::Multiscale::Problem::Six::Diffusion::JacobianRangeType &flux) const {

  // coeff.first = 1.01 + cos( 2.0 * M_PI * (x[0] / epsilon) ) + stochastic perturbation
  // coeff.second = 1.01 + cos( 2.0 * M_PI * (x[0] / epsilon) ) + stochastic perturbation
  const auto coeff = constants().coefficients(x);

  flux[0][0] = coeff.first * direction_gradient[0][0];
  flux[0][1] = coeff.second * direction_gradient[0][1];

}

void Dune::Multiscale::Problem::Six::ExactSolution::evaluate(const Dune::Multiscale::Problem::Six::ExactSolution::DomainType &, Dune::Multiscale::Problem::Six::ExactSolution::RangeType &) const {
  DUNE_THROW(Dune::NotImplemented, "Exact solution not available!");
}

void Dune::Multiscale::Problem::Six::ExactSolution::evaluateJacobian(const Dune::Multiscale::Problem::Six::ExactSolution::DomainType &, Dune::Multiscale::Problem::Six::ExactSolution::JacobianRangeType &) const {
  DUNE_THROW(Dune::NotImplemented, "Exact solution not available!");
}

void Dune::Multiscale::Problem::Six::ExactSolution::evaluate(const Dune::Multiscale::Problem::Six::ExactSolution::DomainType &x, const Dune::Multiscale::Problem::Six::ExactSolution::TimeType &, Dune::Multiscale::Problem::Six::ExactSolution::RangeType &y) const {
  evaluate(x, y);
}
