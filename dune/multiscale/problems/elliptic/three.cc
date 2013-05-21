#include "three.hh"

namespace Dune {
namespace Multiscale {
namespace Problem {
namespace Three {

ModelProblemData::ModelProblemData()
  : IModelProblemData(constants()){
    assert( constants_.epsilon != 0.0);
    if (constants().get("linear", true))
       DUNE_THROW(Dune::InvalidStateException, "problem three is entirely nonlinear, but problem.linear was true");
    if (constants().get("stochastic_pertubation", false) && !(this->problemAllowsStochastics()) )
       DUNE_THROW(Dune::InvalidStateException, "The problem does not allow stochastic perturbations. Please, switch the key off.");
}

inline std::string ModelProblemData::getMacroGridFile() const {
  return("../dune/multiscale/grids/macro_grids/elliptic/earth.dgf");
}

inline bool ModelProblemData::problemIsPeriodic() const {
  return false; // = problem is not periodic
}

inline bool ModelProblemData::problemAllowsStochastics() const {
  return false; // = problem does not allow stochastic perturbations
  // (if you want it, you must add the 'perturb' method provided
  // by 'constants.hh' - see model problems 4 to 7 for examples )
}

} //namespace Three
} //namespace Problem
} //namespace Multiscale {
} //namespace Dune {

