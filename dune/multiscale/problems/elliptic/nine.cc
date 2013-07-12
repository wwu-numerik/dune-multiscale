#include "nine.hh"

namespace Dune {
namespace Multiscale {
namespace Problem {
namespace Nine {

// default value for epsilon (if not sprecified in the parameter file)
CONSTANTSFUNCTION( 0.05 )

ModelProblemData::ModelProblemData()
  : IModelProblemData(constants()) {
    if (!constants().get("linear", true))
      DUNE_THROW(Dune::InvalidStateException, "problem nine is entirely linear, but problem.linear was false");
    if (constants().get("stochastic_pertubation", false) && !(this->problemAllowsStochastics()) )
       DUNE_THROW(Dune::InvalidStateException, "The problem does not allow stochastic perturbations. Please, switch the key off.");
}

std::string ModelProblemData::getMacroGridFile() const {
  return("../dune/multiscale/grids/macro_grids/elliptic/cube_three.dgf");
}

bool ModelProblemData::problemIsPeriodic() const {
  return true; // = problem is periodic
}

bool ModelProblemData::problemAllowsStochastics() const {
  return false; // = problem does not allow stochastic perturbations
  // (if you want it, you must add the 'perturb' method provided
  // by 'constants.hh' - see model problems 4 to 7 for examples )
}

FirstSource::FirstSource(){}

// evaluate f, i.e. return y=f(x) for a given x
// the following method defines 'f':
void __attribute__((hot)) FirstSource::evaluate(const DomainType& x,
                     RangeType& y) const {

  const double eps = constants().epsilon;
  const double pi_square = pow(M_PI, 2.0);
  const double x0_eps = (x[0] / eps);
  const double cos_2_pi_x0_eps = cos( 2.0 * M_PI * x0_eps );
  const double sin_2_pi_x0_eps = sin( 2.0 * M_PI * x0_eps );
  const double coefficient_0 = 2.0 * ( 1.0 / (8.0 * M_PI * M_PI) ) * ( 1.0 / ( 2.0 + cos_2_pi_x0_eps ) );
  const double coefficient_1 = ( 1.0 / (8.0 * M_PI * M_PI) ) * ( 1.0 + ( 0.5 * cos_2_pi_x0_eps ) );
  const double sin_2_pi_x0 = sin(2.0 * M_PI * x[0]);
  const double cos_2_pi_x0 = cos(2.0 * M_PI * x[0]);
  const double sin_2_pi_x1 = sin(2.0 * M_PI * x[1]);

  const double d_x0_coefficient_0
    = pow(2.0 + cos_2_pi_x0_eps, -2.0) * ( 1.0 / (2.0 * M_PI) ) * (1.0 / eps) * sin_2_pi_x0_eps;

  const RangeType grad_u = (2.0* M_PI* cos_2_pi_x0 * sin_2_pi_x1)
                           +( (-1.0) * eps * M_PI * ( sin_2_pi_x0 * sin_2_pi_x1 * sin_2_pi_x0_eps ) )
                           +( M_PI * ( cos_2_pi_x0 * sin_2_pi_x1 * cos_2_pi_x0_eps ) );


  const RangeType d_x0_x0_u =
                    - (4.0 * pi_square * sin_2_pi_x0 * sin_2_pi_x1)
                    - (2.0 * pi_square * ( eps + (1.0 / eps) ) * cos_2_pi_x0 * sin_2_pi_x1 * sin_2_pi_x0_eps)
                    - (4.0 * pi_square * sin_2_pi_x0 * sin_2_pi_x1 * cos_2_pi_x0_eps);

  const RangeType d_x1_x1_u =
                    - (4.0 * pi_square * sin_2_pi_x0 * sin_2_pi_x1)
                    - (2.0 * pi_square * eps * cos_2_pi_x0 * sin_2_pi_x1 * sin_2_pi_x0_eps);

  y = -(d_x0_coefficient_0 * grad_u)
      -(coefficient_0 * d_x0_x0_u)
      -(coefficient_1 * d_x1_x1_u);
} // evaluate

void FirstSource::evaluate(const DomainType& x, const TimeType& /*time*/, RangeType& y) const {
  evaluate(x, y);
}

Diffusion::Diffusion(){}

void Diffusion::diffusiveFlux(const DomainType& x,
                   const JacobianRangeType& direction,
                   JacobianRangeType& flux) const {
  double coefficient_0 = 2.0 * ( 1.0 / (8.0 * M_PI * M_PI) ) * ( 1.0 / ( 2.0 + cos( 2.0 * M_PI * (x[0] / constants().epsilon) ) ) );
  double coefficient_1 = ( 1.0 / (8.0 * M_PI * M_PI) ) * ( 1.0 + ( 0.5 * cos( 2.0 * M_PI * (x[0] / constants().epsilon) ) ) );

  double stab = 0.0;

  flux[0][0] = coefficient_0 * direction[0][0] + stab * direction[0][1];
  flux[0][1] = coefficient_1 * direction[0][1] + stab * direction[0][0];
} // diffusiveFlux

void Diffusion::jacobianDiffusiveFlux(const DomainType& x,
                           const JacobianRangeType& /*position_gradient*/,
                           const JacobianRangeType& direction_gradient,
                           JacobianRangeType& flux) const {
  double coefficient_0 = 2.0 * ( 1.0 / (8.0 * M_PI * M_PI) ) * ( 1.0 / ( 2.0 + cos( 2.0 * M_PI * (x[0] / constants().epsilon) ) ) );
  double coefficient_1 = ( 1.0 / (8.0 * M_PI * M_PI) ) * ( 1.0 + ( 0.5 * cos( 2.0 * M_PI * (x[0] / constants().epsilon) ) ) );
  flux[0][0] = coefficient_0 * direction_gradient[0][0];
  flux[0][1] = coefficient_1 * direction_gradient[0][1];
} // jacobianDiffusiveFlux

ExactSolution::ExactSolution(){}

void ExactSolution::evaluate(const DomainType& x,
                     RangeType& y) const {
  // approximation obtained by homogenized solution + first corrector

  // coarse part
  y = sin(2.0 * M_PI * x[0]) * sin(2.0 * M_PI * x[1]);

  // fine part // || u_fine_part ||_L2 = 0.00883883 (for eps = 0.05 )
  y += 0.5 * constants().epsilon * ( cos(2.0 * M_PI * x[0]) * sin(2.0 * M_PI * x[1]) * sin( 2.0 * M_PI * (x[0] / constants().epsilon) ) );
} // evaluate

void ExactSolution::jacobian(const DomainType& x, typename FunctionSpaceType::JacobianRangeType& grad_u) const {
  grad_u[0][0] = 2.0* M_PI* cos(2.0 * M_PI * x[0]) * sin(2.0 * M_PI * x[1]);
  grad_u[0][1] = 2.0* M_PI* sin(2.0 * M_PI * x[0]) * cos(2.0 * M_PI * x[1]);

  grad_u[0][0] += (-1.0) * constants().epsilon * M_PI
                  * ( sin(2.0 * M_PI * x[0]) * sin(2.0 * M_PI * x[1]) * sin( 2.0 * M_PI * (x[0] / constants().epsilon) ) );
  grad_u[0][0] += M_PI * ( cos(2.0 * M_PI * x[0]) * sin(2.0 * M_PI * x[1]) * cos( 2.0 * M_PI * (x[0] / constants().epsilon) ) );

  grad_u[0][1] += constants().epsilon * M_PI
                  * ( cos(2.0 * M_PI * x[0]) * cos(2.0 * M_PI * x[1]) * sin( 2.0 * M_PI * (x[0] / constants().epsilon) ) );
} // jacobian

void ExactSolution::evaluate(const DomainType& x,
                     const TimeType& /*timedummy*/,
                     RangeType& y) const {
  evaluate(x, y);
}
} //namespace Nine
} //namespace Problem
} //namespace Multiscale {
} //namespace Dune {
