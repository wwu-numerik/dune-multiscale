#include <config.h>
#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/parallel/collectivecommunication.hh>
#include <dune/stuff/common/validation.hh>
#include <dune/stuff/common/float_cmp.hh>
#include <dune/stuff/grid/boundaryinfo.hh>
#include <dune/stuff/common/profiler.hh>
#include <dune/stuff/common/configuration.hh>
#include <dune/multiscale/problems/selector.hh>
#include <math.h>
#include <sstream>
#include <mpi.h>

#include "dune/multiscale/problems/base.hh"

#include "random.hh"
#include "random_permeability.hh"

namespace Dune {
namespace Multiscale {
namespace Problem {
namespace Random {

ModelProblemData::ModelProblemData()
  : IModelProblemData()
  , boundaryInfo_(DSG::BoundaryInfos::NormalBased<typename View::Intersection>::create(boundary_settings()))
  , subBoundaryInfo_() {}

// functions actually unused
std::string ModelProblemData::getMacroGridFile() const {
  return "";
}

std::pair<CommonTraits::DomainType, CommonTraits::DomainType> ModelProblemData::gridCorners() const {
  CommonTraits::DomainType lowerLeft(0.0);
  CommonTraits::DomainType upperRight(1.0);
  return {lowerLeft, upperRight};
}

void ModelProblemData::problem_init(MPIHelper::MPICommunicator global, MPIHelper::MPICommunicator local)
{
  getMutableDiffusion().init(global, local);
}

void Diffusion::init(MPIHelper::MPICommunicator global, MPIHelper::MPICommunicator local)
{
  const int log2Seg = DSC_CONFIG_GET("grids.macro_cells_per_dim", 4);
  int seed = 0;
  MPI_Comm_rank(global, &seed);
  assert(seed >= 0);
  const int overlap = DSC_CONFIG_GET("grids.overlap", 1u);
  correlation_ = DSC::make_unique<CorrelationType>();
  DSC::ScopedTiming field_tm("msfem.perm_field");
  field_ = DSC::make_unique<PermeabilityType>(local, *correlation_, log2Seg, seed+1, overlap);
  field_->create();
}

const ModelProblemData::BoundaryInfoType& ModelProblemData::boundaryInfo() const { return *boundaryInfo_; }

const ModelProblemData::SubBoundaryInfoType& ModelProblemData::subBoundaryInfo() const { return subBoundaryInfo_; }

ParameterTree ModelProblemData::boundary_settings() const {
  Dune::ParameterTree boundarySettings;
  if (DSC_CONFIG.has_sub("problem.boundaryInfo")) {
    boundarySettings = DSC_CONFIG.sub("problem.boundaryInfo");
  } else {
    boundarySettings["default"] = "dirichlet";
    boundarySettings["compare_tolerance"] = "1e-10";
    switch (CommonTraits::world_dim) {
      case 1:
      case 3:
        DUNE_THROW(InvalidStateException, "no boundary settings for 3D random field problem yet");
      case 2:
        boundarySettings["neumann.0"] = "[0.0 1.0]";
        boundarySettings["neumann.1"] = "[0.0 -1.0]";
        break;
    }
  }
  return boundarySettings;
}

Source::Source() {}
Diffusion::Diffusion() {}

PURE HOT void Source::evaluate(const DomainType& /*x*/, RangeType& y) const {
  y = 0.;
}

void Diffusion::evaluate(const DomainType& x, Diffusion::RangeType& ret) const {
   assert(field_);
  const double scalar = field_->operator()(x);
  ret[0][0] = 1 * scalar;
  ret[1][1] = ret[0][0];
}

PURE HOT void Diffusion::diffusiveFlux(const DomainType& x, const Problem::JacobianRangeType& direction,
                                       Problem::JacobianRangeType& flux) const {
  Diffusion::RangeType eval;
  evaluate(x, eval);
  flux[0][0] = eval[0][0] * direction[0][0];
  flux[0][1] = eval[1][1] * direction[0][1];
} // diffusiveFlux

size_t Diffusion::order() const { return 2; }
size_t Source::order() const { return 1; } // evaluate

PURE void DirichletData::evaluate(const DomainType& x, RangeType& y) const {
  y = DSC::FloatCmp::eq(x[0], 1.0) ? 0.0 : 1.0;
} // evaluate

PURE void DirichletData::jacobian(const DomainType& /*x*/, JacobianRangeType& y) const {
  y = 0.;
} // jacobian

PURE void NeumannData::evaluate(const DomainType& /*x*/, RangeType& y) const {
  y = 0.;
} // evaluate

} // namespace Nine
} // namespace Problem
} // namespace Multiscale {
} // namespace Dune {
