#include "eleven.hh"

namespace Dune {
namespace Multiscale {
namespace Problem {
namespace Eleven {

ModelProblemData::ModelProblemData()
  : IModelProblemData(constants()) {
    if (!constants().get("linear", true))
      DUNE_THROW(Dune::InvalidStateException, "problem eleven is entirely linear, but problem.linear was false");
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


// evaluate f, i.e. return y=f(x) for a given x
// the following method defines 'f':
void FirstSource::evaluate(const DomainType& x,
                     RangeType& y) const {

  double coefficient_0 = 2.0 * ( 1.0 / (8.0 * M_PI * M_PI) ) * ( 1.0 / ( 2.0 + cos( 2.0 * M_PI * (x[0] / constants().epsilon) ) ) );
  double coefficient_1 = ( 1.0 / (8.0 * M_PI * M_PI) ) * ( 1.0 + ( 0.5 * cos( 2.0 * M_PI * (x[0] / constants().epsilon) ) ) );

  double d_x0_coefficient_0
    = pow(2.0 + cos( 2.0 * M_PI * (x[0] / constants().epsilon) ), -2.0) * ( 1.0 / (2.0 * M_PI) ) * (1.0 / constants().epsilon) * sin(
    2.0 * M_PI * (x[0] / constants().epsilon) );

  JacobianRangeType grad_u;
  grad_u[0][0] = (20.0) * 2.0 * (1.0 - x[1]) * x[1] * (2.0 * x[0] - 1.0);
  grad_u[0][1] = (20.0) * 2.0 * (1.0 - x[0]) * x[0] * (2.0 * x[1] - 1.0);

  grad_u[0][0] += (-1.0) * constants().epsilon * M_PI
                  * ( sin(2.0 * M_PI * x[0]) * sin(2.0 * M_PI * x[1]) * sin( 2.0 * M_PI * (x[0] / constants().epsilon) ) );
  grad_u[0][0] += M_PI * ( cos(2.0 * M_PI * x[0]) * sin(2.0 * M_PI * x[1]) * cos( 2.0 * M_PI * (x[0] / constants().epsilon) ) );

  grad_u[0][1] += constants().epsilon * M_PI
                  * ( cos(2.0 * M_PI * x[0]) * cos(2.0 * M_PI * x[1]) * sin( 2.0 * M_PI * (x[0] / constants().epsilon) ) );

  RangeType d_x0_x0_u(0.0);
  d_x0_x0_u += (20.0) * 4.0 * (1.0 - x[1]) * x[1];
  d_x0_x0_u -= 2.0
               * pow(M_PI,
                     2.0) * ( constants().epsilon + (1.0 / constants().epsilon) ) * cos(2.0 * M_PI * x[0]) * sin(2.0 * M_PI * x[1]) * sin(
    2.0 * M_PI * (x[0] / constants().epsilon) );
  d_x0_x0_u -= 4.0
               * pow(M_PI, 2.0) * sin(2.0 * M_PI * x[0]) * sin(2.0 * M_PI * x[1]) * cos( 2.0 * M_PI * (x[0] / constants().epsilon) );

  RangeType d_x1_x1_u(0.0);
  d_x1_x1_u += (20.0) * 4.0 * (1.0 - x[0]) * x[0];
  d_x1_x1_u -= 2.0
               * pow(M_PI,
                     2.0) * constants().epsilon
               * cos(2.0 * M_PI * x[0]) * sin(2.0 * M_PI * x[1]) * sin( 2.0 * M_PI * (x[0] / constants().epsilon) );

  y = 0.0;
  y -= d_x0_coefficient_0 * grad_u[0][0];
  y -= coefficient_0 * d_x0_x0_u;
  y -= coefficient_1 * d_x1_x1_u;
  
#if 1
  RangeType u = (20.0) * (-2.0) * (1.0 - x[0]) * (1.0 - x[1]) * x[0] * x[1];
  u += 0.5 * constants().epsilon * ( cos(2.0 * M_PI * x[0]) * sin(2.0 * M_PI * x[1]) * sin( 2.0 * M_PI * (x[0] / constants().epsilon) ) );

  // F(x,u(x), \grad u(x)
  RangeType F_value = 0.0;

#if 1
  double factor = scaling_factor_ * ( 1.0 / (8.0 * M_PI * M_PI) ) * ( 2.0 + cos( 2.0 * M_PI * (x[0] / pow(constants().epsilon, 1.5 ) ) ) );

  double p = u;
  double g_p = 0.0;
    
  double q_1 = 0.25 * sqrt( 8.0 / 7.0 );
  double q_2 = sqrt( 7.0 / 8.0 );  

  // coefficients of the smoothing polynomial q(x) = a x^3 + b x^2 + c x + d
  double a = ( ( 16.0 * q_1 ) + (128.0 * q_2 ) );
  double b = (- 2.0) * ( q_1 - (27.0 / 16.0 ) * a ) ;
  double c = (2.0 * b) - (3.0 * a);
  double d = - (2.0 * a) + b;
    
  if ( p <= -3.0 )
  {
    std::cout << "In FirstSource::evaluate: Nonlinearity not defined for p <= -3, but p = " << p << "." << std::endl;
    abort();
  }
    
  if ( p <= -1.25 )
  { g_p = sqrt( 0.5 * p + 1.5); }
  else if ( p <= -1.0 )
       { g_p = (a * pow( p , 3.0 )) + (b * pow( p , 2.0 )) + (c * p ) + d; }
       else
       { g_p = 0.0; }
  F_value = factor * g_p * grad_u[0][1];
#endif
#if 0
   double scale = 0.1;
   if ( u <= (-0.25) )
   { F_value = -1.0; }
   else if ( u >= 0.25 )
        { F_value = 1.0; }
        else
	{ F_value = sin( 2.0 * M_PI * u); }
   F_value *= (3.0 + ( scale * ( cos( grad_u[0][0] ) + cos( grad_u[0][1] ) ) ) );
#endif
  y += F_value;
#endif
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


// dummy
void Nonlinearity::evaluate(const DomainType& x, RangeType& y) const
{ DUNE_THROW(Dune::NotImplemented, "Inadmissible call for 'evaluate'"); }

// dummy
void Nonlinearity::evaluate(const DomainType& x, const TimeType& /*time*/, RangeType& y) const
{ DUNE_THROW(Dune::NotImplemented, "Inadmissible call for 'evaluate'"); }

// evaluate y = F(x, u(x), \grad u(x))
// 'position = u(x)', 'direction_gradient = \grad u(x)'
void Nonlinearity::evaluate(const DomainType& x,
                            const RangeType& position,
                            const JacobianRangeType& direction_gradient,
                            RangeType& y) const {

   // F(x,p, (z_1,z_2) ) = c g(p) z_2
#if 1

    double scaling_factor = scaling_factor_;
    scaling_factor *= ( 1.0 / (8.0 * M_PI * M_PI) ) * ( 2.0 + cos( 2.0 * M_PI * (x[0] / pow(constants().epsilon, 1.5) ) ) );

    double p = position;
    double g_p = 0.0;
    
    double q_1 = 0.25 * sqrt( 8.0 / 7.0 );
    double q_2 = sqrt( 7.0 / 8.0 );  
  
    // coefficients of the smoothing polynomial q(x) = a x^3 + b x^2 + c x + d
    double a = ( ( 16.0 * q_1 ) + (128.0 * q_2 ) );
    double b = (- 2.0) * ( q_1 - (27.0 / 16.0 ) * a ) ;
    double c = (2.0 * b) - (3.0 * a);
    double d = - (2.0 * a) + b;
    
    if ( p <= -3.0 )
    {
      //!std::cout << "In Nonlinearity::evaluate: Nonlinearity not defined for p <= -3, but p = " << p << "." << std::endl;
      p = -2.999; //!abort();
    }

    if ( p <= -1.25 )
    { g_p = sqrt( 0.5 * p + 1.5); }
    else if ( p <= -1.0 )
         { g_p = (a * pow( p , 3.0 )) + (b * pow( p , 2.0 )) + (c * p ) + d; }
         else
         { g_p = 0.0; }
    y = scaling_factor * g_p * direction_gradient[0][1];
#endif
#if 0
   double scale = 0.1;
   if ( position <= (-0.25) )
   { y = -1.0; }
   else if ( position >= 0.25 )
        { y = 1.0; }
        else
	{ y = sin( 2.0 * M_PI * position); }
   y *= (3.0 + ( scale * ( cos( direction_gradient[0][0] ) + cos( direction_gradient[0][1] ) ) ) );
#endif
}  // evaluate

// evaluate position derivative y = d_1 F (x, u(x), \grad u(x))  (derivative with respect to the second componenent 'u(x)')
// 'position = u(x)', 'direction_gradient = \grad u(x)'
void Nonlinearity::position_derivative(const DomainType& x,
                                       const RangeType& position,
                                       const JacobianRangeType& direction_gradient,
                                       RangeType& y) const {

   // \partial_p F(x,p, (z_1,z_2) ) = ..
#if 1
    double scaling_factor = scaling_factor_;
    scaling_factor *= ( 1.0 / (8.0 * M_PI * M_PI) ) * ( 2.0 + cos( 2.0 * M_PI * (x[0] / pow(constants().epsilon, 1.5) ) ) );

    double p = position;
    double g_p = 0.0;
    
    double q_1 = 0.25 * sqrt( 8.0 / 7.0 );
    double q_2 = sqrt( 7.0 / 8.0 );  
  
    // coefficients of the smoothing polynomial q(x) = a x^3 + b x^2 + c x + d
    double a = ( ( 16.0 * q_1 ) + (128.0 * q_2 ) );
    double b = (- 2.0) * ( q_1 - (27.0 / 16.0 ) * a ) ;
    double c = (2.0 * b) - (3.0 * a);

    if ( p <= -3.0 )
    {
      //!std::cout << "In Nonlinearity::position_derivative: Nonlinearity not defined for p <= -3, but p = " << p << "." << std::endl;
      p = -2.999; //!abort();
    }
 
    if ( p <= -1.25 )
    { g_p = 0.25 * pow( 0.5 * p + 1.5, -0.5 ); }
    else if ( p <= -1.0 )
         { g_p = (3.0 * a * pow( p , 2.0 )) + (2.0 * b * p) + c; }
         else
         { g_p = 0.0; }
    y = scaling_factor * g_p * direction_gradient[0][1];
#endif
#if 0
   double scale = 0.1;
   if ( position <= (-0.25) )
   { y = 0.0; }
   else if ( position >= 0.25 )
        { y = 0.0; }
        else
	{ y = 2.0 * M_PI * cos( 2.0 * M_PI * position); }
   y *= (3.0 + ( scale * ( cos( direction_gradient[0][0] ) + cos( direction_gradient[0][1] ) ) ) );
#endif
}  // position_derivative

// evaluate position derivative y = d_2 F (x, u(x), \grad u(x))  (derivative with respect to the third componenent 'grad u(x)')
// 'position = u(x)', 'direction_gradient = \grad u(x)'
void Nonlinearity::direction_derivative(const DomainType& x,
                                        const RangeType& position,
                                        const JacobianRangeType& direction_gradient,
                                        JacobianRangeType& y) const {
   // \grad_z F(x,p, (z_1,z_2) ) = ...
#if 1
    double scaling_factor = scaling_factor_;
    scaling_factor *= ( 1.0 / (8.0 * M_PI * M_PI) ) * ( 2.0 + cos( 2.0 * M_PI * (x[0] / pow(constants().epsilon, 1.5) ) ) );

    double p = position;
    double g_p = 0.0;
    
    double q_1 = 0.25 * sqrt( 8.0 / 7.0 );
    double q_2 = sqrt( 7.0 / 8.0 );  
  
    // coefficients of the smoothing polynomial q(x) = a x^3 + b x^2 + c x + d
    double a = ( ( 16.0 * q_1 ) + (128.0 * q_2 ) );
    double b = (- 2.0) * ( q_1 - (27.0 / 16.0 ) * a ) ;
    double c = (2.0 * b) - (3.0 * a);
    double d = - (2.0 * a) + b;
    
    if ( p <= -3.0 )
    {
      //!std::cout << "Nonlinearity::direction_derivative: Nonlinearity not defined for p <= -3. But p = " << p << "." << std::endl;
      p = -2.999; //!abort();
    }
    
    if ( p <= -1.25 )
    { g_p = sqrt( 0.5 * p + 1.5); }
    else if ( p <= -1.0 )
         { g_p = (a * pow( p , 3.0 )) + (b * pow( p , 2.0 )) + (c * p ) + d; }
         else
         { g_p = 0.0; }
    y[0][1] = scaling_factor * g_p;
    y[0][0] = 0.0;
#endif
#if 0
   double scale = 0.1;
   if ( position <= (-0.25) )
   { y[0][0] = -1.0; }
   else if ( position >= 0.25 )
        { y[0][0] = 1.0; }
        else
        { y[0][0] = sin( 2.0 * M_PI * position); }
   y[0][1] = y[0][0];
   y[0][0] *= scale * (-1.0) * sin( direction_gradient[0][0] );
   y[0][1] *= scale * (-1.0) * sin( direction_gradient[0][1] );
#endif

}  // direction_derivative

ExactSolution::ExactSolution(){}

void ExactSolution::evaluate(const DomainType& x,
                     RangeType& y) const {
  // approximation obtained by homogenized solution + first corrector

  // coarse part
  y = (20.0) * (-2.0) * (1.0 - x[0]) * (1.0 - x[1]) * x[0] * x[1]; //sin(2.0 * M_PI * x[0]) * sin(2.0 * M_PI * x[1]);

  // fine part // || u_fine_part ||_L2 = 0.00883883 (for eps = 0.05 )
  y += 0.5 * constants().epsilon * ( cos(2.0 * M_PI * x[0]) * sin(2.0 * M_PI * x[1]) * sin( 2.0 * M_PI * (x[0] / constants().epsilon) ) );
} // evaluate

void ExactSolution::evaluateJacobian(const DomainType& x, typename FunctionSpaceType::JacobianRangeType& grad_u) const {
  grad_u[0][0] = (20.0) * 2.0 * (1.0 - x[1]) * x[1] * (2.0 * x[0] - 1.0);
  grad_u[0][1] = (20.0) * 2.0 * (1.0 - x[0]) * x[0] * (2.0 * x[1] - 1.0);

  grad_u[0][0] += (-1.0) * constants().epsilon * M_PI
                  * ( sin(2.0 * M_PI * x[0]) * sin(2.0 * M_PI * x[1]) * sin( 2.0 * M_PI * (x[0] / constants().epsilon) ) );
  grad_u[0][0] += M_PI * ( cos(2.0 * M_PI * x[0]) * sin(2.0 * M_PI * x[1]) * cos( 2.0 * M_PI * (x[0] / constants().epsilon) ) );

  grad_u[0][1] += constants().epsilon * M_PI
                  * ( cos(2.0 * M_PI * x[0]) * cos(2.0 * M_PI * x[1]) * sin( 2.0 * M_PI * (x[0] / constants().epsilon) ) );
} // evaluateJacobian

void ExactSolution::evaluate(const DomainType& x,
                     const TimeType& /*timedummy*/,
                     RangeType& y) const {
  evaluate(x, y);
}
} //namespace Eleven
} //namespace Problem
} //namespace Multiscale {
} //namespace Dune {
