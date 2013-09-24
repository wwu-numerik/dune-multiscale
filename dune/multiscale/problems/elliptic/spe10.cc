#include "spe10.hh"

#include <boost/assert.hpp>

namespace Dune {
namespace Multiscale {
namespace Problem {
namespace SPE10 {

// default value for epsilon (if not specified in the parameter file)
CONSTANTSFUNCTION( 0.05 )

ModelProblemData::ModelProblemData()
  : IModelProblemData(constants()) {
    if (!constants().get("linear", true))
      DUNE_THROW(Dune::InvalidStateException, "problem SPE10 is entirely linear, but problem.linear was false");
    if (constants().get("stochastic_pertubation", false) && !(this->problemAllowsStochastics()) )
       DUNE_THROW(Dune::InvalidStateException, "The problem does not allow stochastic perturbations. Please, switch the key off.");
}

std::string ModelProblemData::getMacroGridFile() const {
  return("../dune/multiscale/grids/macro_grids/elliptic/spe10.dgf");
}

bool ModelProblemData::problemIsPeriodic() const {
  return false; // = problem is periodic
}

bool ModelProblemData::problemAllowsStochastics() const {
  return false; // = problem does not allow stochastic perturbations
  // (if you want it, you must add the 'perturb' method provided
  // by 'constants.hh' - see model problems 4 to 7 for examples )
}

FirstSource::FirstSource(){}

// evaluate f, i.e. return y=f(x) for a given x
// the following method defines 'f':
void __attribute__((hot)) FirstSource::evaluate(const DomainType& /*x*/,
        RangeType& y) const {
  y = RangeType(1.0);
} // evaluate

void FirstSource::evaluate(const DomainType& x, const TimeType& /*time*/, RangeType& y) const {
  evaluate(x, y);
}


Diffusion::Diffusion() :deltas_{6.096, 3.048, 0.6096},
permeability_(nullptr)
{
  readPermeability();
}

Diffusion::~Diffusion() {
  delete permeability_;
  permeability_=nullptr;
}

void Diffusion::diffusiveFlux(const DomainType& x,
                   const JacobianRangeType& direction,
                   JacobianRangeType& flux) const {
  assert(x.size()<=3 && "SPE 10 model is only defined for up to three dimensions!");
  Dune::FieldMatrix< double, DomainType::dimension, DomainType::dimension> permeability(0.0);

  //TODO this class does not seem to work in 2D, when changing 'spe10.dgf' to a 2D grid?
  
  // 3 is the maximum space dimension
  std::vector<double> intervalls(3, 0.0);
  for (int dim=0; dim<DomainType::dimension; ++dim)
    intervalls[dim] = std::floor(x[dim] / deltas_[dim]);

  int offset = 0;
  switch ( DomainType::dimension ) {
    case 1:
      offset = intervalls[0];
      permeability[0][0] = permeability_[offset];
      break;
    case 2:
      offset = intervalls[0] + intervalls[1] * 60;
      permeability[0][0] = permeability_[offset];
      permeability[1][1] = permeability_[offset + 1122000];
      break;
    case 3:
      offset = intervalls[0] + intervalls[1] * 60 + intervalls[2] * 220 * 60;
      permeability[0][0] = permeability_[offset];
      permeability[1][1] = permeability_[offset + 1122000];
      permeability[2][2] = permeability_[offset + 2244000];
      break;
  }

  permeability.mv(direction[0], flux[0]);
} // diffusiveFlux

void Diffusion::jacobianDiffusiveFlux(const DomainType& /*x*/,
                           const JacobianRangeType& /*position_gradient*/,
                           const JacobianRangeType& /*direction_gradient*/,
                           JacobianRangeType& /*flux*/) const {
  DUNE_THROW(NotImplemented, "Jacobian of Flux is not implemented at the moment!");
} // jacobianDiffusiveFlux

void Diffusion::readPermeability()   {
  permeability_ = new double[3366000];
  std::string   filename = "../dune/multiscale/problems/elliptic/spe10_permeability.dat";
  std::ifstream file(filename.c_str());
  double        val;
  if (!file) { // file couldn't be opened
    DSC_LOG_ERROR << "The SPE10-permeability data file could not be opened. This file does\n"
            << "not come with the dune-multiscale repository due to file size. To download it\n"
            << "execute\n"
            << "wget http://www.spe.org/web/csp/datasets/por_perm_case2a.zip\n"
            << "unzip the file and move the file 'spe_perm.dat' to\n"
            << "dune-multiscale/dune/multiscale/problems/elliptic/spe10_permeability.dat!\n";
    DUNE_THROW(IOError, "Data file for Groundwaterflow permeability could not be opened!");
  }
  file >> val;
  int counter = 0;
  while (!file.eof()) {
    // keep reading until end-of-file
    permeability_[counter++] = val;
    file >> val;         // sets EOF flag if no value found
  }
  file.close();
  return;
}         /* readPermeability */

} //namespace SPE10
} //namespace Problem
} //namespace Multiscale {
} //namespace Dune {
