#include "twelve.hh"

namespace Dune {
namespace Multiscale {
namespace Problem {
namespace Twelve {

// default value for epsilon (if not sprecified in the parameter file)
CONSTANTSFUNCTION( 0.05 )

ModelProblemData::ModelProblemData()
  : IModelProblemData(constants()) {
    if (!constants().get("linear", true))
      DUNE_THROW(Dune::InvalidStateException, "problem Twelve is entirely linear, but problem.linear was false");
    if (constants().get("stochastic_pertubation", false) && !(this->problemAllowsStochastics()) )
       DUNE_THROW(Dune::InvalidStateException, "The problem does not allow stochastic perturbations. Please, switch the key off.");
}

std::string ModelProblemData::getMacroGridFile() const {
  return("../dune/multiscale/grids/macro_grids/elliptic/cube_three.dgf");
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
void FirstSource::evaluate(const DomainType& x,
                     RangeType& y) const
{

  double center = 0.5; // i.e. 'center = (center,center)'
  double radius = 0.05;

  if ( sqrt( pow(x[0]-center,2.0) + pow(x[1]-center,2.0) ) <= radius )
    y = 20.0;
  else
    y = 0.0;
} // evaluate

void FirstSource::evaluate(const DomainType& x, const TimeType& /*time*/, RangeType& y) const {
  evaluate(x, y);
}

Diffusion::Diffusion(){}

void Diffusion::diffusiveFlux(const DomainType& x,
                   const JacobianRangeType& direction,
                   JacobianRangeType& flux) const {

  double isolater_thickness = 0.02;
  double isolater_conductivity = 1e-2;
  double coefficient = 0.0;
  
  bool x_0_qualifies = false;
  bool x_1_qualifies = false;

  if ( (x[0] >= isolater_thickness) && (x[0] <= (2.0*isolater_thickness) ) )
    x_0_qualifies = true;

  if ( ( x[0] >= (1.0-2.0*isolater_thickness) ) && (x[0] <= (1.0 - isolater_thickness) ) )
    x_0_qualifies = true;

  if ( (x[1] >= isolater_thickness) && (x[1] <= (2.0*isolater_thickness) ) )
    x_1_qualifies = true;

  if ( ( x[1] >= (1.0-2.0*isolater_thickness) ) && (x[1] <= (1.0 - isolater_thickness) ) )
    x_1_qualifies = true;
  
  if ( x_0_qualifies )
  {
    if ( (x[1] >= isolater_thickness) && (x[1] <= (1.0 - isolater_thickness) ) )
     x_1_qualifies = true;
    else
     x_1_qualifies = false;
  }

  if ( x_1_qualifies )
  {
    if ( (x[0] >= isolater_thickness) && (x[0] <= (1.0 - isolater_thickness) ) )
     x_0_qualifies = true;
    else
     x_0_qualifies = false;
  }
  
  double eps_0 = constants().epsilon;
  if ( x_0_qualifies && x_1_qualifies )
    coefficient = isolater_conductivity;
  else
  {
      const double x0_eps = (x[0] / eps_0);
      const double cos_2_pi_x0_eps = cos( 2.0 * M_PI * x0_eps );
      coefficient = 0.1 * ( 2.0 + cos_2_pi_x0_eps );
  }
  
  double center = 0.5; // i.e. 'center = (center,center)'
  double radius = 0.25;

  bool high_conductivty = true;
  double eps = eps_0;
  if ( sqrt( pow(x[0]-center,2.0) + pow(x[1]-center,2.0) ) <= radius )
  {
    coefficient = 1.0;
    while ( (radius-eps) >= eps_0 )
    {
      high_conductivty = !high_conductivty;
      if ( sqrt( pow(x[0]-center,2.0) + pow(x[1]-center,2.0) ) <= (radius-eps) )
      {
        if ( high_conductivty )
          coefficient = 1.0;
        else
          coefficient = 1e-1;
      }
      eps += eps_0;
    }
  }
  //coefficient *= 0.001;
  
  double a_0_0 = coefficient;
  double a_0_1 = 0.0;
  double a_1_0 = 0.0;
  double a_1_1 = coefficient;
  
  flux[0][0] = a_0_0 * direction[0][0] + a_0_1 * direction[0][1];
  flux[0][1] = a_1_0 * direction[0][0] + a_1_1 * direction[0][1];
} // diffusiveFlux

void Diffusion::jacobianDiffusiveFlux(const DomainType& x,
                           const JacobianRangeType& /*position_gradient*/,
                           const JacobianRangeType& direction_gradient,
                           JacobianRangeType& flux) const {
  double isolater_thickness = 0.02;
  double isolater_conductivity = 1e-1;
  double coefficient = 0.0;
  
  bool x_0_qualifies = false;
  bool x_1_qualifies = false;

  if ( (x[0] >= isolater_thickness) && (x[0] <= (2.0*isolater_thickness) ) )
    x_0_qualifies = true;

  if ( ( x[0] >= (1.0-2.0*isolater_thickness) ) && (x[0] <= (1.0 - isolater_thickness) ) )
    x_0_qualifies = true;

  if ( (x[1] >= isolater_thickness) && (x[1] <= (2.0*isolater_thickness) ) )
    x_1_qualifies = true;

  if ( ( x[1] >= (1.0-2.0*isolater_thickness) ) && (x[1] <= (1.0 - isolater_thickness) ) )
    x_1_qualifies = true;
  
  if ( x_0_qualifies )
  {
    if ( (x[1] >= isolater_thickness) && (x[1] <= (1.0 - isolater_thickness) ) )
     x_1_qualifies = true;
    else
     x_1_qualifies = false;
  }

  if ( x_1_qualifies )
  {
    if ( (x[0] >= isolater_thickness) && (x[0] <= (1.0 - isolater_thickness) ) )
     x_0_qualifies = true;
    else
     x_0_qualifies = false;
  }
  
  double eps_0 = constants().epsilon;
  if ( x_0_qualifies && x_1_qualifies )
    coefficient = isolater_conductivity;
  else
  {
      const double x0_eps = (x[0] / eps_0);
      const double cos_2_pi_x0_eps = cos( 2.0 * M_PI * x0_eps );
      coefficient = 0.1 * ( 2.0 + cos_2_pi_x0_eps );
  }
  
  double center = 0.5; // i.e. 'center = (center,center)'
  double radius = 0.25;

  bool high_conductivty = true;
  double eps = eps_0;
  if ( sqrt( pow(x[0]-center,2.0) + pow(x[1]-center,2.0) ) <= radius )
  {
    coefficient = 1.0;
    while ( (radius-eps) >= eps_0 )
    {
      high_conductivty = !high_conductivty;
      if ( sqrt( pow(x[0]-center,2.0) + pow(x[1]-center,2.0) ) <= (radius-eps) )
      {
        if ( high_conductivty )
          coefficient = 1.0;
        else
          coefficient = 1e-1;
      }
      eps += eps_0;
    }
  }
  //coefficient *= 0.001;

  double a_0_0 = coefficient;
  double a_0_1 = 0.0;
  double a_1_0 = 0.0;
  double a_1_1 = coefficient;
 
  flux[0][0] = a_0_0 * direction_gradient[0][0] + a_0_1 * direction_gradient[0][1];
  flux[0][1] = a_1_0 * direction_gradient[0][0] + a_1_1 * direction_gradient[0][1];
} // jacobianDiffusiveFlux


// evaluate Dirichlet Boundary Function
void DirichletBoundaryCondition::evaluate(const DomainType& x,
                                          RangeType& y) const
{
  y = x[0];

} // evaluate

void DirichletBoundaryCondition::evaluate(const DomainType& x, const TimeType& /*time*/, RangeType& y) const {
  evaluate(x, y);
}


void Dune::Multiscale::Problem::Twelve::ExactSolution::evaluate(const Dune::Multiscale::Problem::Twelve::ExactSolution::DomainType &, Dune::Multiscale::Problem::Twelve::ExactSolution::RangeType &) const {
  DUNE_THROW(Dune::NotImplemented, "Exact solution not available!");
}

void Dune::Multiscale::Problem::Twelve::ExactSolution::jacobian(const Dune::Multiscale::Problem::Twelve::ExactSolution::DomainType &, Dune::Multiscale::Problem::Twelve::ExactSolution::JacobianRangeType &) const {
  DUNE_THROW(Dune::NotImplemented, "Exact solution not available!");
}

void Dune::Multiscale::Problem::Twelve::ExactSolution::evaluate(const Dune::Multiscale::Problem::Twelve::ExactSolution::DomainType &x, const Dune::Multiscale::Problem::Twelve::ExactSolution::TimeType &, Dune::Multiscale::Problem::Twelve::ExactSolution::RangeType &y) const {
  evaluate(x, y);
}

} //namespace Twelve
} //namespace Problem
} //namespace Multiscale {
} //namespace Dune {
