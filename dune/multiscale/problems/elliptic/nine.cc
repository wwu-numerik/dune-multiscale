#include "nine.hh"

namespace Dune {
namespace Multiscale {
namespace Problem {
namespace Nine {

ModelProblemData::ModelProblemData()
  : IModelProblemData(constants()) {
    if (!constants().get("linear", true))
      DUNE_THROW(Dune::InvalidStateException, "problem nine is entirely linear, but problem.linear was false");
    if (constants().get("stochastic_pertubation", false) && !(this->problemAllowsStochastics()) )
       DUNE_THROW(Dune::InvalidStateException, "The problem does not allow stochastic perturbations. Please, switch the key off.");
}

inline std::string ModelProblemData::getMacroGridFile() const {
  return("../dune/multiscale/grids/macro_grids/elliptic/msfem_cube_three.dgf");
}

inline bool ModelProblemData::problemIsPeriodic() const {
  return true; // = problem is periodic
}

inline bool ModelProblemData::problemAllowsStochastics() const {
  return false; // = problem does not allow stochastic perturbations
  // (if you want it, you must add the 'perturb' method provided
  // by 'constants.hh' - see model problems 4 to 7 for examples )
}

} //namespace Nine
} //namespace Problem
} //namespace Multiscale {
} //namespace Dune {
