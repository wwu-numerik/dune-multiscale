#include "one.hh"

namespace Dune {
namespace Multiscale {
namespace Problem {
namespace One {

ModelProblemData::ModelProblemData()
  : IModelProblemData(constants()){
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
  return false; // = problem does not allow stochastic perturbations
  // (if you want it, you must add the 'perturb' method provided
  // by 'constants.hh' - see model problems 4 to 7 for examples )
}

} //namespace One
} //namespace Problem
} //namespace Multiscale {
} //namespace Dune {

