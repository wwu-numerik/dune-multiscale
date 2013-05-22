#include "five.hh"

namespace Dune {
namespace Multiscale {
namespace Problem {
namespace Five {

ModelProblemData::ModelProblemData()
  : IModelProblemData(constants()) {
    assert( constants_.epsilon != 0.0);
    if (constants().get("linear", true))
      DUNE_THROW(Dune::InvalidStateException, "Problem five is entirely nonlinear, but problem.linear was true.");
    if (constants().get("stochastic_pertubation", false) && !(this->problemAllowsStochastics()) )
       DUNE_THROW(Dune::InvalidStateException, "The problem does not allow stochastic perturbations. Please, switch the key off.");
}

std::string ModelProblemData::getMacroGridFile() const {
  return("../dune/multiscale/grids/macro_grids/elliptic/corner_singularity.dgf");
}

bool ModelProblemData::problemIsPeriodic() const {
  return true; // = problem is periodic
}

bool ModelProblemData::problemAllowsStochastics() const {
  return true; // = problem allows stochastic perturbations
}

} //namespace Five
} //namespace Problem
} //namespace Multiscale {
} //namespace Dune {

