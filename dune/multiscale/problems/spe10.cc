#include <config.h>
#include <assert.h>
#include <dune/common/exceptions.hh>
#include <dune/xt/common/logging.hh>
#include <dune/xt/common/validation.hh>
#include <dune/xt/common/configuration.hh>
#include <sstream>

#include "dune/multiscale/problems/base.hh"
#include "spe10.hh"

namespace Dune {
namespace Multiscale {
namespace Problem {
namespace SPE10 {

static const std::string model2_filename = "spe10_permeability.dat";
static const size_t model2_x_elements = 60;
static const size_t model2_y_elements = 220;
static const size_t model2_z_elements = 85;
static const double model_2_length_x = 365.76;
static const double model_2_length_y = 670.56;
static const double model_2_length_z = 51.816;
static const double model_2_min_value = 6.65e-8; // isotropic: 0.000665
static const double model_2_max_value = 20000;

ModelProblemData::ModelProblemData(MPIHelper::MPICommunicator global,
                                   MPIHelper::MPICommunicator local,
                                   Dune::XT::Common::Configuration config_in)
  : IModelProblemData(global, local, config_in)
  , subBoundaryInfo_()
{
  boundaryInfo_ = std::unique_ptr<ModelProblemData::BoundaryInfoType>(
      Dune::XT::Grid::BoundaryInfos::NormalBased<typename View::Intersection>::create(boundary_settings()));
  subBoundaryInfo_ = std::unique_ptr<ModelProblemData::SubBoundaryInfoType>(
      Dune::XT::Grid::BoundaryInfos::NormalBased<typename SubView::Intersection>::create(boundary_settings()));
}

Dune::XT::Common::Configuration default_config()
{
  Dune::XT::Common::Configuration config;
  config["filename"] = model2_filename;
  config["name"] = "Spe10Diffusion";
  config["lower_left"] = "[0 0 0]";
  config["upper_right"] = "[" + Dune::XT::Common::to_string(model_2_length_x) + " "
                          + Dune::XT::Common::to_string(model_2_length_y) + " "
                          + Dune::XT::Common::to_string(model_2_length_z) + "]";
  config["upper_right"] = "[2 5 1]";
  config["anisotropic"] = "true";
  config["min"] = Dune::XT::Common::to_string(model_2_min_value);
  config["max"] = Dune::XT::Common::to_string(model_2_max_value);
  return config;
} // ... default_config(...)

std::pair<CommonTraits::DomainType, CommonTraits::DomainType> ModelProblemData::gridCorners() const
{
  CommonTraits::DomainType lowerLeft(0.0);
  CommonTraits::DomainType upperRight(0.0);
  if (View::dimension != 3)
    DUNE_THROW(Dune::InvalidStateException, "SPE data only available for world dim == 3");
  const auto ll = default_config().template get<CommonTraits::DomainType>("lower_left");
  const auto ur = default_config().template get<CommonTraits::DomainType>("upper_right");
  return {ll, ur};
}

const ModelProblemData::BoundaryInfoType& ModelProblemData::boundaryInfo() const
{
  return *boundaryInfo_;
}

const ModelProblemData::SubBoundaryInfoType& ModelProblemData::subBoundaryInfo() const
{
  return *subBoundaryInfo_;
}

ParameterTree ModelProblemData::boundary_settings() const
{
  Dune::ParameterTree boundarySettings;
  if (DXTC_CONFIG.has_sub("problem.boundaryInfo")) {
    boundarySettings = DXTC_CONFIG.sub("problem.boundaryInfo");
  } else {
    boundarySettings["default"] = "neumann";
    boundarySettings["compare_tolerance"] = "1e-10";
    switch (View::dimension /*View is defined in IModelProblemData*/) {
      case 1:
        DUNE_THROW(NotImplemented, "Boundary values are not implemented for SPE10 in 1D!");
        break;
      case 2:
        boundarySettings["dirichlet.0"] = "[0.0 1.0]";
        break;
      case 3:
        boundarySettings["dirichlet.0"] = "[0.0 1.0 0.0]";
    }
  }
  return boundarySettings;
}

void DirichletData::evaluate(const DomainType& /*x*/, RangeType& y) const
{
  y = 0.0;
} // evaluate

void NeumannData::evaluate(const DomainType& x, RangeType& y) const
{
  if (Dune::XT::Common::FloatCmp::eq(x[1], CommonTraits::RangeFieldType(0)))
    y = 1.0;
  else
    y = 0.0;
} // evaluate

Source::Source(MPIHelper::MPICommunicator /*global*/,
               MPIHelper::MPICommunicator /*local*/,
               Dune::XT::Common::Configuration /*config_in*/)
{
}

void __attribute__((hot)) Source::evaluate(const DomainType& /*x*/, RangeType& y) const
{
  y = typename CommonTraits::RangeType(0.0);
} // evaluate

Diffusion::Diffusion(MPIHelper::MPICommunicator /*global*/,
                     MPIHelper::MPICommunicator /*local*/,
                     Dune::XT::Common::Configuration /*config_in*/)
  : permeability_(nullptr)
{
  readPermeability();
}

Diffusion::~Diffusion()
{
}

void Diffusion::evaluate(const DomainType& x, Diffusion::RangeType& y) const
{
  BOOST_ASSERT_MSG(x.size() == 3, "SPE 10 model is only defined for three dimensions!");
  // TODO this class does not seem to work in 2D, when changing 'spe10.dgf' to a 2D grid?
  if (!permeability_) {
    MS_LOG_ERROR_0 << "The SPE10-permeability data file could not be opened. This file does\n"
                   << "not come with the dune-multiscale repository due to file size. To download it\n"
                   << "execute\n"
                   << "wget http://www.spe.org/web/csp/datasets/por_perm_case2a.zip\n"
                   << "unzip the file and move the file 'spe_perm.dat' to\n"
                   << "dune-multiscale/dune/multiscale/problems/spe10_permeability.dat!\n";
    DUNE_THROW(IOError, "Data file for Groundwaterflow permeability could not be opened!");
  }
  // 3 is the maximum space dimension

  const auto center = x;
  std::vector<size_t> whichPartition(dimDomain, 0);
  const Dune::XT::Common::FieldVector<CommonTraits::DomainFieldType, CommonTraits::world_dim> ll =
      default_config().template get<CommonTraits::DomainType>("lower_left");
  const Dune::XT::Common::FieldVector<CommonTraits::DomainFieldType, CommonTraits::world_dim> ur =
      default_config().template get<CommonTraits::DomainType>("upper_right");
  // code need to compile, not run, for dim 2, therefore hardcode to 3 elements
  std::array<size_t, 3> ne{{model2_x_elements, model2_y_elements, model2_z_elements}};

  for (size_t dd = 0; dd < dimDomain; ++dd) {
    // for points that are on upperRight_[d], this selects one partition too much
    // so we need to cap this
    whichPartition[dd] = std::min(size_t(std::floor(ne[dd] * ((center[dd] - ll[dd]) / (ur[dd] - ll[dd])))), ne[dd] - 1);
  }
  size_t subdomain = 0;
  if (dimDomain == 1)
    subdomain = whichPartition[0];
  else if (dimDomain == 2)
    subdomain = whichPartition[0] + whichPartition[1] * ne[0];
  else
    subdomain = whichPartition[0] + whichPartition[1] * ne[0] + whichPartition[2] * ne[1] * ne[0];
  const auto ii = subdomain;
  const size_t entries_per_dim = model2_x_elements * model2_y_elements * model2_z_elements;
  const bool anisotropic = true;
  y[0][0] = permeability_[ii];
  y[1][1] = permeability_[anisotropic ? entries_per_dim + ii : ii];
  y[2][2] = permeability_[anisotropic ? 2 * entries_per_dim + ii : ii];
}

void Diffusion::diffusiveFlux(const DomainType& x,
                              const Problem::JacobianRangeType& direction,
                              Problem::JacobianRangeType& flux) const
{
  Diffusion::RangeType eval_tmp;
  evaluate(x, eval_tmp);
  eval_tmp.mv(direction[0], flux[0]);
} // diffusiveFlux

void Diffusion::readPermeability()
{
  auto min = default_config().template get<RangeFieldType>("min");
  auto max = default_config().template get<RangeFieldType>("max");
  std::ifstream datafile(model2_filename);
  if (!datafile.is_open())
    DUNE_THROW(spe10_model2_data_file_missing, "could not open '" << model2_filename << "'!");
  if (!(max > min))
    DUNE_THROW(Dune::RangeError, "max (is " << max << ") has to be larger than min (is " << min << ")!");
  const RangeFieldType scale = (max - min) / (model_2_max_value - model_2_min_value);
  const RangeFieldType shift = min - scale * model_2_min_value;
  // read all the data from the file
  const size_t entries_per_dim = model2_x_elements * model2_y_elements * model2_z_elements;
  //  std::vector< double > data(3*entries_per_dim);
  const size_t data_size = 3 * entries_per_dim;
  permeability_ = std::unique_ptr<double[]>(new double[data_size]());
  double tmp = 0;
  size_t counter = 0;
  while (datafile >> tmp) {
    permeability_[counter++] = (tmp * scale) + shift;
  }
  datafile.close();
  if (counter != data_size)
    DUNE_THROW(Dune::IOError,
               "wrong number of entries in '" << model2_filename << "' (are " << counter << ", should be " << data_size
                                              << ")!");
  return;
} /* readPermeability */

} // namespace SPE10
} // namespace Problem
} // namespace Multiscale {
} // namespace Dune {
