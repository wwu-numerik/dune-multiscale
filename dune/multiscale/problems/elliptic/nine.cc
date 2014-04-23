#include <config.h>
#include <dune/common/exceptions.hh>
#include <dune/stuff/common/parameter/validation.hh>
#include <dune/stuff/grid/boundaryinfo.hh>
#include <math.h>
#include <sstream>

#include "dune/multiscale/problems/base.hh"
#include "dune/multiscale/problems/constants.hh"
#include "nine.hh"

namespace Dune {
namespace Multiscale {
namespace Problem {
namespace Nine {

// default value for epsilon (if not specified in the parameter file)
CONSTANTSFUNCTION(0.05)

ModelProblemData::ModelProblemData()
  : IModelProblemData(constants())
  , boundaryInfo_(Stuff::GridboundaryNormalBased<typename View::Intersection>::create(boundary_settings()))
  , subBoundaryInfo_()
{}

std::string ModelProblemData::getMacroGridFile() const {
  return ("../dune/multiscale/grids/macro_grids/elliptic/msfem_cube_three.dgf");
}

bool ModelProblemData::problemIsPeriodic() const { return true; }

bool ModelProblemData::problemAllowsStochastics() const { return false; }

std::pair<CommonTraits::DomainType, CommonTraits::DomainType>
ModelProblemData::gridCorners() const {
  CommonTraits::DomainType lowerLeft(0.0);
  CommonTraits::DomainType upperRight(1.0);
  return {lowerLeft, upperRight};
}

const ModelProblemData::BoundaryInfoType& ModelProblemData::boundaryInfo() const {
  return *boundaryInfo_;
}

const ModelProblemData::SubBoundaryInfoType& ModelProblemData::subBoundaryInfo() const {
  return subBoundaryInfo_;
}

ParameterTree ModelProblemData::boundary_settings() const {
  Dune::ParameterTree boundarySettings;
  if (DSC_CONFIG.hasSub("problem.boundaryInfo")) {
    boundarySettings = DSC_CONFIG.sub("problem.boundaryInfo");
  } else {
    boundarySettings["default"] = "dirichlet";
    boundarySettings["compare_tolerance"] = "1e-10";
    switch (View::dimension /*View is defined in IModelProblemData*/) {
    case 1:
      break;
    case 2:
      break;
    case 3:
      boundarySettings["neumann.0"] = "[0.0; 0.0; 1.0]";
      boundarySettings["neumann.1"] = "[0.0; 0.0; -1.0]";
    }
  }
  return boundarySettings;
}


FirstSource::FirstSource() {}

void __attribute__((hot)) FirstSource::evaluate(const DomainType& x, RangeType& y) const {

  const double eps = constants().epsilon;
  const double pi_square = pow(M_PI, 2.0);
  const double x0_eps = (x[0] / eps);
  const double cos_2_pi_x0_eps = cos(2.0 * M_PI * x0_eps);
  const double sin_2_pi_x0_eps = sin(2.0 * M_PI * x0_eps);
  const double coefficient_0 = 2.0 * (1.0 / (8.0 * M_PI * M_PI)) * (1.0 / (2.0 + cos_2_pi_x0_eps));
  const double coefficient_1 = (1.0 / (8.0 * M_PI * M_PI)) * (1.0 + (0.5 * cos_2_pi_x0_eps));
  const double sin_2_pi_x0 = sin(2.0 * M_PI * x[0]);
  const double cos_2_pi_x0 = cos(2.0 * M_PI * x[0]);
  const double sin_2_pi_x1 = sin(2.0 * M_PI * x[1]);

  const double d_x0_coefficient_0 =
      pow(2.0 + cos_2_pi_x0_eps, -2.0) * (1.0 / (2.0 * M_PI)) * (1.0 / eps) * sin_2_pi_x0_eps;

  const typename FunctionSpaceType::RangeType grad_u =
      (2.0 * M_PI * cos_2_pi_x0 * sin_2_pi_x1) + ((-1.0) * eps * M_PI * (sin_2_pi_x0 * sin_2_pi_x1 * sin_2_pi_x0_eps)) +
      (M_PI * (cos_2_pi_x0 * sin_2_pi_x1 * cos_2_pi_x0_eps));

  const typename FunctionSpaceType::RangeType d_x0_x0_u =
      -(4.0 * pi_square * sin_2_pi_x0 * sin_2_pi_x1) -
      (2.0 * pi_square * (eps + (1.0 / eps)) * cos_2_pi_x0 * sin_2_pi_x1 * sin_2_pi_x0_eps) -
      (4.0 * pi_square * sin_2_pi_x0 * sin_2_pi_x1 * cos_2_pi_x0_eps);

  const typename FunctionSpaceType::RangeType d_x1_x1_u =
      -(4.0 * pi_square * sin_2_pi_x0 * sin_2_pi_x1) -
      (2.0 * pi_square * eps * cos_2_pi_x0 * sin_2_pi_x1 * sin_2_pi_x0_eps);

  y = -(d_x0_coefficient_0 * grad_u) - (coefficient_0 * d_x0_x0_u) - (coefficient_1 * d_x1_x1_u);
} // evaluate

void FirstSource::evaluate(const DomainType& x, const TimeType& /*time*/, RangeType& y) const { evaluate(x, y); }

Diffusion::Diffusion() {}

void Diffusion::diffusiveFlux(const DomainType& x, const JacobianRangeType& direction, JacobianRangeType& flux) const {
  double coefficient_0 =
      2.0 * (1.0 / (8.0 * M_PI * M_PI)) * (1.0 / (2.0 + cos(2.0 * M_PI * (x[0] / constants().epsilon))));
  double coefficient_1 = (1.0 / (8.0 * M_PI * M_PI)) * (1.0 + (0.5 * cos(2.0 * M_PI * (x[0] / constants().epsilon))));

  double stab = 0.0;

  flux[0][0] = coefficient_0 * direction[0][0] + stab * direction[0][1];
  flux[0][1] = coefficient_1 * direction[0][1] + stab * direction[0][0];
} // diffusiveFlux

void Diffusion::jacobianDiffusiveFlux(const DomainType& x, const JacobianRangeType& /*position_gradient*/,
                                      const JacobianRangeType& direction_gradient, JacobianRangeType& flux) const {
  double coefficient_0 =
      2.0 * (1.0 / (8.0 * M_PI * M_PI)) * (1.0 / (2.0 + cos(2.0 * M_PI * (x[0] / constants().epsilon))));
  double coefficient_1 = (1.0 / (8.0 * M_PI * M_PI)) * (1.0 + (0.5 * cos(2.0 * M_PI * (x[0] / constants().epsilon))));
  flux[0][0] = coefficient_0 * direction_gradient[0][0];
  flux[0][1] = coefficient_1 * direction_gradient[0][1];
} // jacobianDiffusiveFlux

ExactSolution::ExactSolution() {}

void ExactSolution::evaluate(const DomainType& x, RangeType& y) const {
  // approximation obtained by homogenized solution + first corrector

  const double eps = constants().epsilon;
  static const double M_TWOPI = M_PI * 2.0;
  const double x0_eps = (x[0] / eps);
  const double sin_2_pi_x0_eps = sin(M_TWOPI * x0_eps);
  const double x0_2_pi = M_TWOPI * x[0];
  const double x1_2_pi = M_TWOPI * x[1];
  const double sin_2_pi_x0 = sin(x0_2_pi);
  const double cos_2_pi_x0 = cos(x0_2_pi);
  const double sin_2_pi_x1 = sin(x1_2_pi);

  y = sin_2_pi_x0 * sin_2_pi_x1
      + (0.5 * eps * cos_2_pi_x0 * sin_2_pi_x1 * sin_2_pi_x0_eps);
} // evaluate

void __attribute__((hot)) ExactSolution::jacobian(const DomainType& x, JacobianRangeType& grad_u) const {
  const double eps = constants().epsilon;
  static const double M_TWOPI = M_PI * 2.0;
  const double x0_eps = (x[0] / eps);
  const double cos_2_pi_x0_eps = cos(M_TWOPI * x0_eps);
  const double sin_2_pi_x0_eps = sin(M_TWOPI * x0_eps);
  const double x0_2_pi = M_TWOPI * x[0];
  const double x1_2_pi = M_TWOPI * x[1];
  const double sin_2_pi_x0 = sin(x0_2_pi);
  const double cos_2_pi_x0 = cos(x0_2_pi);
  const double sin_2_pi_x1 = sin(x1_2_pi);
  const double cos_2_pi_x1 = cos(x1_2_pi);

  grad_u[0][1] = (M_TWOPI * sin_2_pi_x0 * cos_2_pi_x1) + (eps * M_PI * cos_2_pi_x0 * cos_2_pi_x1 * sin_2_pi_x0_eps);

  grad_u[0][0] = (M_TWOPI * cos_2_pi_x0 * sin_2_pi_x1) - (eps * M_PI * sin_2_pi_x0 * sin_2_pi_x1 * sin_2_pi_x0_eps) +
                 (M_PI * cos_2_pi_x0 * sin_2_pi_x1 * cos_2_pi_x0_eps);

} // jacobian

void DirichletData::evaluate(const DomainType& /*x*/, RangeType& y) const { y = 0.0; } // evaluate

void DirichletData::evaluate(const DomainType& x, const TimeType& /*time*/, RangeType& y) const { evaluate(x, y); }

void DirichletData::jacobian(const DomainType& /*x*/, JacobianRangeType& y) const { y[0] = 0.0; } // jacobian

void DirichletData::jacobian(const DomainType& x, const TimeType& /*time*/, JacobianRangeType& y) const {
  jacobian(x, y);
} // jacobian

void NeumannData::evaluate(const DomainType& /*x*/, RangeType& y) const { y = 0.0; } // evaluate

void NeumannData::evaluate(const DomainType& x, const TimeType& /*time*/, RangeType& y) const { evaluate(x, y); }

} // namespace Nine
} // namespace Problem
} // namespace Multiscale {
} // namespace Dune {
