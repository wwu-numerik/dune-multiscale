#include "ten.hh"

namespace Dune {
namespace Multiscale {
namespace Problem {
namespace Ten {

ModelProblemData::ModelProblemData()
  : IModelProblemData(constants()) {
  assert(constants_.epsilon == 0.0);
  if (!constants().get("linear", true))
    DUNE_THROW(Dune::InvalidStateException, "problem ten is entirely linear, but problem.linear was false");
  if (constants().get("stochastic_pertubation", false) && !(this->problemAllowsStochastics()) )
     DUNE_THROW(Dune::InvalidStateException, "The problem does not allow stochastic perturbations. Please, switch the key off.");

}

std::string ModelProblemData::getMacroGridFile() const {
  return("../dune/multiscale/grids/macro_grids/elliptic/cube_three.dgf");      // _strange_grid
}

bool ModelProblemData::problemIsPeriodic() const {
  return false; // = problem is not periodic
}

bool ModelProblemData::problemAllowsStochastics() const {
  return false; // = problem does not allow stochastic perturbations
  // (if you want it, you must add the 'perturb' method provided
  // by 'constants.hh' - see model problems 4 to 7 for examples )
}

} //namespace Ten
} //namespace Problem
} //namespace Multiscale {
} //namespace Dune {

