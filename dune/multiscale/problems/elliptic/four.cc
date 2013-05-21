#include "four.hh"

namespace Dune {
namespace Multiscale {
namespace Problem {
namespace Four {

ModelProblemData::ModelProblemData()
  : IModelProblemData(constants()) {
    assert( constants_.epsilon != 0.0);
    if (constants().get("stochastic_pertubation", false) && !(this->problemAllowsStochastics()) )
       DUNE_THROW(Dune::InvalidStateException, "The problem does not allow stochastic perturbations. Please, switch the key off.");
}

inline std::string ModelProblemData::getMacroGridFile() const {
  return("../dune/multiscale/grids/macro_grids/elliptic/corner_singularity.dgf");
}

inline bool ModelProblemData::problemIsPeriodic() const {
  return true; // = problem is periodic
}

inline bool ModelProblemData::problemAllowsStochastics() const {
  return true; // = problem allows stochastic perturbations
}

} //namespace Four
} //namespace Problem
} //namespace Multiscale {
} //namespace Dune {
