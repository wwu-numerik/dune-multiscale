#include "example.hh"

namespace Dune {
namespace Multiscale {
namespace Problem {
namespace Example {

// default value for epsilon (=0.05)
// (in case epsilon is not specified in the paramter file)
CONSTANTSFUNCTION(0.05) // 0.05 is a dummy! no epsilon in our example

ModelProblemData::ModelProblemData() : IModelProblemData(constants()) {
  if (constants().get("stochastic_pertubation", false) && !(this->problemAllowsStochastics()))
    DUNE_THROW(Dune::InvalidStateException,
               "The problem does not allow stochastic perturbations. Please, switch the key off.");
}

std::string ModelProblemData::getMacroGridFile() const {
  return ("../dune/multiscale/grids/macro_grids/elliptic/cube_one.dgf"); // a standard 2d-cube
}

bool ModelProblemData::problemIsPeriodic() const {
  return false; // = problem is not periodic
}

bool ModelProblemData::problemAllowsStochastics() const {
  return false; // = problem does not allow stochastic perturbations
                // (if you want it, you must add the 'perturb' method provided
                // by 'constants.hh' - see model problems 4 to 7 for examples )
}

} // namespace Example
} // namespace Problem
} // namespace Multiscale {
} // namespace Dune {

void Dune::Multiscale::Problem::Example::FirstSource::evaluate(
    const Dune::Multiscale::Problem::Example::FirstSource::DomainType& x,
    Dune::Multiscale::Problem::Example::FirstSource::RangeType& y) const {

  // f(x_1,x_2) :=  sin( 2 π x_1 ) · sin( 2 π x_2 )
  y = sin(2.0 * M_PI * x[0]) * sin(2.0 * M_PI * x[1]);
}

void Dune::Multiscale::Problem::Example::Diffusion::diffusiveFlux(
    const Dune::Multiscale::Problem::Example::Diffusion::DomainType&,
    const Dune::Multiscale::Problem::Example::Diffusion::JacobianRangeType& direction,
    Dune::Multiscale::Problem::Example::Diffusion::JacobianRangeType& flux) const {

  // coefficient = ( 1 /(8π²) )
  double coefficient = 1.0 * (1.0 / (8.0 * M_PI * M_PI));
  //  for 'direction = ∇v(x)', we set 'flux = A(x,∇v(x)) = ( 1 /(8π²) ) ∇v(x)':
  flux[0][0] = coefficient * direction[0][0];
  flux[0][1] = coefficient * direction[0][1];
}

void Dune::Multiscale::Problem::Example::Diffusion::jacobianDiffusiveFlux(
    const Dune::Multiscale::Problem::Example::Diffusion::DomainType&,
    const Dune::Multiscale::Problem::Example::Diffusion::JacobianRangeType&,
    const Dune::Multiscale::Problem::Example::Diffusion::JacobianRangeType& direction_gradient,
    Dune::Multiscale::Problem::Example::Diffusion::JacobianRangeType& flux) const {
  // Note: in the linear case, we have
  //     'JA( x , position_gradient ) = A( x )'
  // for all directions.

  // coefficient = ( 1 /(8π²) )
  double coefficient = 1.0 * (1.0 / (8.0 * M_PI * M_PI));
  flux[0][0] = coefficient * direction_gradient[0][0];
  flux[0][1] = coefficient * direction_gradient[0][1];
}

void Dune::Multiscale::Problem::Example::ExactSolution::evaluate(
    const Dune::Multiscale::Problem::Example::ExactSolution::DomainType& x,
    Dune::Multiscale::Problem::Example::ExactSolution::RangeType& y) const {
  // u( x_1, x_2 ) = sin( 2 π x_1 ) · sin( 2 π x_2 )
  y = sin(2.0 * M_PI * x[0]) * sin(2.0 * M_PI * x[1]);
}

void Dune::Multiscale::Problem::Example::ExactSolution::evaluate(
    const Dune::Multiscale::Problem::Example::ExactSolution::DomainType& x,
    const Dune::Multiscale::Problem::Example::ExactSolution::TimeType&,
    Dune::Multiscale::Problem::Example::ExactSolution::RangeType& y) const {
  evaluate(x, y);
}
