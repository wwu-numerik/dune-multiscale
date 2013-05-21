#include "six.hh"

namespace Dune {
namespace Multiscale {
namespace Problem {
namespace Six {

ModelProblemData::ModelProblemData()
  : IModelProblemData(constants()) {
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


} //namespace Six
} //namespace Problem
} //namespace Multiscale {
} //namespace Dune {
