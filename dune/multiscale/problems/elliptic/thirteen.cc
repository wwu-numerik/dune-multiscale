#include <config.h>
#include <dune/common/exceptions.hh>
#include <dune/stuff/common/parameter/validation.hh>
#include <sstream>

#include "dune/multiscale/problems/base.hh"
#include "dune/multiscale/problems/constants.hh"
#include "thirteen.hh"

namespace Dune {
namespace Multiscale {
namespace Problem {
namespace Thirteen {

// default value for epsilon (if not sprecified in the parameter file)
CONSTANTSFUNCTION(0.05)

ModelProblemData::ModelProblemData() : IModelProblemData(constants()) {
  if (constants().get("stochastic_pertubation", false) && !(this->problemAllowsStochastics()))
    DUNE_THROW(Dune::InvalidStateException,
               "The problem does not allow stochastic perturbations. Please, switch the key off.");
}

std::string ModelProblemData::getMacroGridFile() const {
  return ("../dune/multiscale/grids/macro_grids/elliptic/cube_three_dirichlet_neumann.dgf");
}

bool ModelProblemData::problemIsPeriodic() const { return false; }

bool ModelProblemData::problemAllowsStochastics() const { return false; }

void FirstSource::evaluate(const DomainType& /*x*/, RangeType& y) const { y = 0.0; } // evaluate

void FirstSource::evaluate(const DomainType& x, const TimeType& /*time*/, RangeType& y) const { evaluate(x, y); }

Diffusion::Diffusion() {}

void Diffusion::diffusiveFlux(const DomainType& x, const JacobianRangeType& direction, JacobianRangeType& flux) const {

  double conductor_thickness = 0.05;
  double conductivity = 20.0;

  double isolator_thickness = 0.05;
  double isolator_conductivity = 1e-2;

  double coefficient = 0.0;

  bool in_conductor = false;
  bool in_isolator = false;

  bool x_0_qualifies = false;
  bool x_1_qualifies = false;

  // two stribes of high conductivity
  double position_x_a = 0.0;
  double position_x_b = 0.8;

  double position_y_a = 0.2;
  double position_y_b = 0.8 - conductor_thickness;

  if ((x[0] >= position_x_a) && (x[0] <= position_x_b))
    x_0_qualifies = true;

  if (x_0_qualifies) {
    if ((x[1] >= position_y_a) && (x[1] <= (position_y_a + conductor_thickness))) {
      //      coefficient = conductivity*(x[0]-position_x_a)*(position_x_b-x[0]);
      //      coefficient *= (x[1]-position_y_a)*((position_y_a+conductor_thickness)-x[1]);
      //      coefficient /= conductor_thickness;
      //      coefficient /= conductor_thickness;
      //      coefficient += 1.0;
      coefficient = conductivity;
      in_conductor = true;
    }

    if ((x[1] >= position_y_b) && (x[1] <= (position_y_b + conductor_thickness))) {
      //      coefficient = conductivity*(x[0]-position_x_a)*(position_x_b-x[0]);
      //      coefficient *= (x[1]-position_y_b)*((position_y_b+conductor_thickness)-x[1]);
      //      coefficient /= conductor_thickness;
      //      coefficient /= conductor_thickness;
      //      coefficient += 1.0;
      coefficient = conductivity;
      in_conductor = true;
    }
  }

  x_0_qualifies = false;
  x_1_qualifies = false;

  if ((x[0] >= (0.5 - 0.5 * isolator_thickness)) && (x[0] <= (0.5 + 0.5 * isolator_thickness)))
    x_0_qualifies = true;

  if ((x[1] >= 0.35) && (x[1] <= 0.65))
    x_1_qualifies = true;

  if (x_0_qualifies && x_1_qualifies) {
    coefficient = isolator_conductivity;
    in_isolator = true;
  }

  if (!in_isolator && !in_conductor) {
    double eps = constants().epsilon;
    double floor_x_0 = std::floor(x[0] + x[1]);
    double floor_x_1 = std::floor(x[1] - x[0]);
    coefficient = 1.2 + (0.5 * sin(floor_x_0 + std::floor(x[0] / eps) + std::floor(x[1] / eps)));
    coefficient += (0.5 * cos(floor_x_1 + std::floor(x[0] / eps) + std::floor(x[1] / eps)));
    // coefficient = 1.0;
  }

  double a_0_0 = coefficient;
  double a_0_1 = 0.0;
  double a_1_0 = 0.0;
  double a_1_1 = coefficient;

  flux[0][0] = a_0_0 * direction[0][0] + a_0_1 * direction[0][1];
  flux[0][1] = a_1_0 * direction[0][0] + a_1_1 * direction[0][1];
} // diffusiveFlux

void Diffusion::jacobianDiffusiveFlux(const DomainType& x, const JacobianRangeType& /*position_gradient*/,
                                      const JacobianRangeType& direction_gradient, JacobianRangeType& flux) const {

  double conductor_thickness = 0.05;
  double conductivity = 100.0;

  double isolator_thickness = 0.05;
  double isolator_conductivity = 1e-2;

  double coefficient = 0.0;

  bool in_conductor = false;
  bool in_isolator = false;

  bool x_0_qualifies = false;
  bool x_1_qualifies = false;

  // two stribes of high conductivity
  double position_x_a = 0.0;
  double position_x_b = 0.8;

  double position_y_a = 0.2;
  double position_y_b = 0.8 - conductor_thickness;

  if ((x[0] >= position_x_a) && (x[0] <= position_x_b))
    x_0_qualifies = true;

  if (x_0_qualifies) {
    if ((x[1] >= position_y_a) && (x[1] <= (position_y_a + conductor_thickness))) {
      //      coefficient = conductivity*(x[0]-position_x_a)*(position_x_b-x[0]);
      //      coefficient *= (x[1]-position_y_a)*((position_y_a+conductor_thickness)-x[1]);
      //      coefficient /= conductor_thickness;
      //      coefficient /= conductor_thickness;
      //      coefficient += 1.0;
      coefficient = conductivity;
      in_conductor = true;
    }

    if ((x[1] >= position_y_b) && (x[1] <= (position_y_b + conductor_thickness))) {
      //      coefficient = conductivity*(x[0]-position_x_a)*(position_x_b-x[0]);
      //      coefficient *= (x[1]-position_y_b)*((position_y_b+conductor_thickness)-x[1]);
      //      coefficient /= conductor_thickness;
      //      coefficient /= conductor_thickness;
      //      coefficient += 1.0;
      coefficient = conductivity;
      in_conductor = true;
    }
  }

  x_0_qualifies = false;
  x_1_qualifies = false;

  if ((x[0] >= (0.5 - 0.5 * isolator_thickness)) && (x[0] <= (0.5 + 0.5 * isolator_thickness)))
    x_0_qualifies = true;

  if ((x[1] >= 0.35) && (x[1] <= 0.65))
    x_1_qualifies = true;

  if (x_0_qualifies && x_1_qualifies) {
    coefficient = isolator_conductivity;
    in_isolator = true;
  }

  if (!in_isolator && !in_conductor) {
    double eps = constants().epsilon;
    double floor_x_0 = std::floor(x[0] + x[1]);
    double floor_x_1 = std::floor(x[1] - x[0]);
    coefficient = 1.2 + (0.5 * sin(floor_x_0 + std::floor(x[0] / eps) + std::floor(x[1] / eps)));
    coefficient += (0.5 * cos(floor_x_1 + std::floor(x[0] / eps) + std::floor(x[1] / eps)));
    // coefficient = 1.0;
  }

  double a_0_0 = coefficient;
  double a_0_1 = 0.0;
  double a_1_0 = 0.0;
  double a_1_1 = coefficient;

  flux[0][0] = a_0_0 * direction_gradient[0][0] + a_0_1 * direction_gradient[0][1];
  flux[0][1] = a_1_0 * direction_gradient[0][0] + a_1_1 * direction_gradient[0][1];
} // jacobianDiffusiveFlux


void NeumannBoundaryCondition::evaluate(const DomainType& x, RangeType& y) const {

  double conductor_thickness = 0.05;

  double position_y_a = 0.2;
  double position_y_b = 0.8 - conductor_thickness;

  if ((x[1] >= position_y_a) && (x[1] <= (position_y_a + conductor_thickness))) {
    y = 2.0;
  } else if ((x[1] >= position_y_b) && (x[1] <= (position_y_b + conductor_thickness))) {
    y = 2.0;
  } else
    y = 0.0;

} // evaluate

void NeumannBoundaryCondition::evaluate(const DomainType& x, const TimeType& /*time*/, RangeType& y) const {
  evaluate(x, y);
}

} // namespace Thirteen
} // namespace Problem
} // namespace Multiscale {
} // namespace Dune {
