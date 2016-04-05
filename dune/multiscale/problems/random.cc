#include <config.h>
#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/parallel/collectivecommunication.hh>
#include <dune/xt/common/validation.hh>
#include <dune/xt/common/float_cmp.hh>
#include <dune/stuff/grid/boundaryinfo.hh>
#include <dune/xt/common/timings.hh>
#include <dune/xt/common/configuration.hh>
#include <dune/multiscale/problems/selector.hh>
#include <math.h>
#include <sstream>

#include "dune/multiscale/problems/base.hh"

#include "random.hh"

#if HAVE_RANDOM_PROBLEM
#include <mpi.h>
#include "random_permeability.hh"
#include <mpi.h>
#endif

namespace Dune {
namespace Multiscale {
namespace Problem {
namespace Random {

ModelProblemData::ModelProblemData()
  : IModelProblemData()
  , boundaryInfo_(DSG::BoundaryInfos::NormalBased<typename View::Intersection>::create(boundary_settings()))
  , subBoundaryInfo_() {}

// functions actually unused
std::string ModelProblemData::getMacroGridFile() const { return ""; }

std::pair<CommonTraits::DomainType, CommonTraits::DomainType> ModelProblemData::gridCorners() const {
  CommonTraits::DomainType lowerLeft(0.0);
  CommonTraits::DomainType upperRight(1.0);
  return {lowerLeft, upperRight};
}

void ModelProblemData::problem_init(DMP::ProblemContainer& problem, MPIHelper::MPICommunicator global,
                                    MPIHelper::MPICommunicator local) {
  problem.getMutableDiffusion().init(problem, global, local);
}

void ModelProblemData::prepare_new_evaluation(DMP::ProblemContainer& problem) {
  problem.getMutableDiffusion().prepare_new_evaluation();
}

void Diffusion::init(const DMP::ProblemContainer& problem, MPIHelper::MPICommunicator global,
                     MPIHelper::MPICommunicator local) {
  const auto cells_per_dim = problem.config().get<std::vector<std::size_t>>("grids.macro_cells_per_dim");
  std::for_each(cells_per_dim.begin(), cells_per_dim.end(), [&](size_t t) { assert(t == cells_per_dim[0]); });
  const int log2Seg = std::log2l(cells_per_dim[0]);
  int seed = 0;
#if HAVE_RANDOM_PROBLEM
  MPI_Comm_rank(global, &seed);
  assert(seed >= 0);
  const int overlap = problem.config().get("grids.overlap", 1u);
  const auto corrLen = problem.config().get("mlmc.correlation_length", 0.2f);
  const auto sigma = problem.config().get("mlmc.correlation_sigma", 1.0f);
  correlation_ = Dune::XT::Common::make_unique<Correlation>(corrLen, sigma);
  Dune::XT::Common::ScopedTiming field_tm("msfem.perm_field.init");
  field_ = Dune::XT::Common::make_unique<PermeabilityType>(local, *correlation_, log2Seg, seed + 1, overlap);
#else
  DUNE_THROW(InvalidStateException, "random problem needs additional libs to be configured properly");
#endif
}

void Diffusion::prepare_new_evaluation() {
  Dune::XT::Common::ScopedTiming field_tm("msfem.perm_field.create");
#if HAVE_RANDOM_PROBLEM
  assert(field_);
  field_->create();
#else
  DUNE_THROW(InvalidStateException, "random problem needs additional libs to be configured properly");
#endif
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
        DUNE_THROW(InvalidStateException, "no boundary settings for 1D random field problem");
      case 3:
        boundarySettings["neumann.0"] = "[0.0  1.0  0.0]";
        boundarySettings["neumann.1"] = "[0.0 -1.0  0.0]";
        boundarySettings["neumann.2"] = "[0.0  0.0  1.0]";
        boundarySettings["neumann.3"] = "[0.0  0.0 -1.0]";
        break;
      case 2:
        boundarySettings["neumann.0"] = "[0.0  1.0]";
        boundarySettings["neumann.1"] = "[0.0 -1.0]";
        break;
    }
  }
  return boundarySettings;
}

Diffusion::Diffusion() {}

void Diffusion::evaluate(const DomainType& x, Diffusion::RangeType& ret) const {
#if HAVE_RANDOM_PROBLEM
  assert(field_);
  const double scalar = field_->operator()(x);
  ret = 0;
  for (const auto i : Dune::XT::Common::value_range(CommonTraits::world_dim))
    ret[i][i] = scalar;
#else
  DUNE_THROW(InvalidStateException, "random problem needs additional libs to be configured properly");
#endif
}

PURE HOT void Diffusion::diffusiveFlux(const DomainType& x, const Problem::JacobianRangeType& direction,
                                       Problem::JacobianRangeType& flux) const {
  Diffusion::RangeType eval;
  evaluate(x, eval);
  flux = 0;
  for (const auto i : Dune::XT::Common::value_range(CommonTraits::world_dim))
    flux[0][i] = eval[i][i] * direction[0][i];

} // diffusiveFlux

size_t Diffusion::order() const { return 1; }

PURE void DirichletData::evaluate(const DomainType& x, RangeType& y) const { y = 1.0 - x[0]; } // evaluate

PURE void NeumannData::evaluate(const DomainType& /*x*/, RangeType& y) const { y = 0.; } // evaluate

} // namespace Nine
} // namespace Problem
} // namespace Multiscale {
} // namespace Dune {
