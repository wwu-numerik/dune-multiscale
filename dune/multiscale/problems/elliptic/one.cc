#include <config.h>
#include <assert.h>
#include <dune/common/exceptions.hh>
#include <dune/stuff/common/parameter/validation.hh>
#include <math.h>
#include <sstream>

#include "dune/multiscale/problems/base.hh"
#include "dune/multiscale/problems/constants.hh"
#include "one.hh"

namespace Dune {
namespace Multiscale {
namespace Problem {
namespace One {

// default value for epsilon (if not sprecified in the parameter file)
CONSTANTSFUNCTION(0.05)

ModelProblemData::ModelProblemData() : IModelProblemData(constants()) {
  assert(constants_.epsilon != 0.0);
  if (constants().get("stochastic_pertubation", false) && !(this->problemAllowsStochastics()))
    DUNE_THROW(Dune::InvalidStateException,
               "The problem does not allow stochastic perturbations. Please, switch the key off.");
}

std::string ModelProblemData::getMacroGridFile() const {
  return ("../dune/multiscale/grids/macro_grids/elliptic/cube_two.dgf");
}

bool ModelProblemData::problemIsPeriodic() const {
  return true; // = problem is periodic
}

bool ModelProblemData::problemAllowsStochastics() const { return false; }

void Diffusion::diffusiveFlux(const DomainType& x, const JacobianRangeType& gradient, JacobianRangeType& flux) const {

  const auto diffusion_coefficient = 2.0 + sin(2.0 * M_PI * (x[0] / constants().epsilon));

  flux[0][0] = diffusion_coefficient * gradient[0][0];
  flux[0][1] = diffusion_coefficient * gradient[0][1];
}

void Diffusion::jacobianDiffusiveFlux(const DomainType& x, const JacobianRangeType&,
                                      const JacobianRangeType& direction_gradient, JacobianRangeType& flux) const {

  const auto diffusion_coefficient = 2.0 + sin(2.0 * M_PI * (x[0] / constants().epsilon));

  flux[0][0] = diffusion_coefficient * direction_gradient[0][0];
  flux[0][1] = diffusion_coefficient * direction_gradient[0][1];
}

} // namespace One
} // namespace Problem
} // namespace Multiscale {
} // namespace Dune {
