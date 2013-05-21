#include "seven.hh"

namespace Dune {
namespace Multiscale {
namespace Problem {
namespace Seven {

ModelProblemData::ModelProblemData()
  : IModelProblemData(constants()) {
    if (constants().get("linear", true))
      DUNE_THROW(Dune::InvalidStateException, "Problem seven is entirely nonlinear, but problem.linear was true.");
    assert( constants_.epsilon != 0.0);
    if (constants().get("stochastic_pertubation", false) && !(this->problemAllowsStochastics()) )
       DUNE_THROW(Dune::InvalidStateException, "The problem does not allow stochastic perturbations. Please, switch the key off.");
  }

inline std::string ModelProblemData::getMacroGridFile() const {
  return("../dune/multiscale/grids/macro_grids/elliptic/cube_two.dgf");
}

inline bool ModelProblemData::problemIsPeriodic() const {
  return true; // = problem is periodic
}

inline bool ModelProblemData::problemAllowsStochastics() const {
  return true; // = problem allows stochastic perturbations
}


} //namespace Seven
} //namespace Problem
} //namespace Multiscale {
} //namespace Dune {

