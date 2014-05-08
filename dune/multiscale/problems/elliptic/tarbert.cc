#include <config.h>
#include <assert.h>
#include <dune/common/exceptions.hh>
#include <dune/stuff/common/logging.hh>
#include <dune/stuff/common/parameter/validation.hh>
#include <dune/stuff/common/parameter/configcontainer.hh>
#include <sstream>

#include "dune/multiscale/problems/base.hh"
#include "tarbert.hh"


namespace Dune {
namespace Multiscale {
namespace Problem {
namespace Tarbert {


ModelProblemData::ModelProblemData()
  : IModelProblemData()
  , subBoundaryInfo_()
{
  boundaryInfo_ = std::unique_ptr<ModelProblemData::BoundaryInfoType>(DSG::BoundaryInfos::NormalBased<typename View::Intersection>::create(boundary_settings()));
  subBoundaryInfo_ = std::unique_ptr<ModelProblemData::SubBoundaryInfoType>(DSG::BoundaryInfos::NormalBased<typename SubView::Intersection>::create(boundary_settings()));
}

std::string ModelProblemData::getMacroGridFile() const {
  return ("../dune/multiscale/grids/macro_grids/elliptic/spe10.dgf");
}

bool ModelProblemData::problemIsPeriodic() const { return false; }

bool ModelProblemData::problemAllowsStochastics() const { return false; }

std::pair<CommonTraits::DomainType, CommonTraits::DomainType>
ModelProblemData::gridCorners() const {
  CommonTraits::DomainType lowerLeft(0.0);
  CommonTraits::DomainType upperRight(0.0);
  switch (View::dimension /*View is defined in IModelProblemData*/) {
    case 1:
      DUNE_THROW(NotImplemented, "Tarbert is not defined in 1D!");
      break;
    case 2:
      upperRight[0] = 365.76;
      upperRight[1] = 670.56;
      break;
    case 3:
      upperRight[0] = 365.76;
      upperRight[1] = 670.56;
      upperRight[2] = 21.336;
  }
  return {lowerLeft, upperRight};
}

const ModelProblemData::BoundaryInfoType& ModelProblemData::boundaryInfo() const {
  return *boundaryInfo_;
}

const ModelProblemData::SubBoundaryInfoType& ModelProblemData::subBoundaryInfo() const {
  return *subBoundaryInfo_;
}

ParameterTree ModelProblemData::boundary_settings() const {
  Dune::ParameterTree boundarySettings;
  if (DSC_CONFIG.hasSub("problem.boundaryInfo")) {
    boundarySettings = DSC_CONFIG.sub("problem.boundaryInfo");
  } else {
    boundarySettings["default"] = "neumann";
    boundarySettings["compare_tolerance"] = "1e-10";
    switch (View::dimension /*View is defined in IModelProblemData*/) {
      case 1:
        DUNE_THROW(NotImplemented, "Boundary values are not implemented for Tarbert in 1D!");
        break;
      case 2:
        boundarySettings["dirichlet.0"] = "[0.0; -1.0]";
        break;
      case 3:
        boundarySettings["dirichlet.0"] = "[0.0; 1.0; 0.0]";
//        boundarySettings["dirichlet.1"] = "[0.0; 1.0; 0.0]";
    }
  }
  return boundarySettings;
}

void DirichletData::evaluate(const DomainType& /*x*/, RangeType& y) const {
  y = 0.0;
} // evaluate

void NeumannData::evaluate(const DomainType& x, RangeType& y) const {
  if (std::abs(x[1])<1e-6)
    y = -1.0;
  else
    y = 0.0;
} // evaluate

Source::Source() {}

void __attribute__((hot)) Source::evaluate(const DomainType& /*x*/, RangeType& y) const {
  y = typename FunctionSpaceType::RangeType(0.0);
} // evaluate

Diffusion::Diffusion()
  : deltas_{6.096, 3.048, 0.6096}
  , permeability_(nullptr)
  , permMatrix_(0.0) {
  readPermeability();
}

Diffusion::~Diffusion() {
  delete permeability_;
  permeability_ = nullptr;
}

void Diffusion::diffusiveFlux(const DomainType& x, const JacobianRangeType& direction, JacobianRangeType& flux) const {
  BOOST_ASSERT_MSG(x.size() <= 3, "Tarbert model is only defined for up to three dimensions!");
  // TODO this class does not seem to work in 2D, when changing 'spe10.dgf' to a 2D grid?

  // 3 is the maximum space dimension
  for (int dim = 0; dim < DomainType::dimension; ++dim)
    permIntervalls_[dim] = std::floor(x[dim] / deltas_[dim]);

  int offset = 0;
  switch (DomainType::dimension) {
    case 1:
      DUNE_THROW(NotImplemented, "Tarbert is not implemented for 1D!");
    case 2:
      offset = permIntervalls_[0] + permIntervalls_[1] * 60;
      permMatrix_[0][0] = permeability_[offset];
      permMatrix_[1][1] = permeability_[offset + 1122000];
      break;
    case 3:
      // for the tarbert formation, we always skip the first 50 layers
      offset = permIntervalls_[0] + permIntervalls_[1] * 60 + (permIntervalls_[2]) * 220 * 60;
      permMatrix_[0][0] = permeability_[offset];
      permMatrix_[1][1] = permeability_[offset + 1122000];
      permMatrix_[2][2] = permeability_[offset + 2244000];
      break;
  }

  permMatrix_.mv(direction[0], flux[0]);
//  flux[0]=direction[0];
} // diffusiveFlux

void Diffusion::jacobianDiffusiveFlux(const DomainType& /*x*/, const JacobianRangeType& /*position_gradient*/,
                                      const JacobianRangeType& /*direction_gradient*/,
                                      JacobianRangeType& /*flux*/) const {
  DUNE_THROW(NotImplemented, "Jacobian of Flux is not implemented at the moment!");
} // jacobianDiffusiveFlux

void Diffusion::readPermeability() {
  permeability_ = new double[3366000];
  std::string filename = "../dune/multiscale/problems/elliptic/spe10_permeability.dat";
  std::ifstream file(filename.c_str());
  double val;
  if (!file) { // file couldn't be opened
    DSC_LOG_INFO_0 << "The SPE10-permeability data file could not be opened. This file does\n"
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
    file >> val; // sets EOF flag if no value found
  }
  file.close();
  return;
} /* readPermeability */
  
//  void Diffusion::visualizePermeability(const CommonTraits::GridType& grid) const
//  {
//    typedef typename CommonTraits::GridType::LeafGridView::template Codim< 0 >::Iterator IteratorType;
//    typedef typename CommonTraits::GridType::LeafGridView::template Codim< 0 >::Entity   EntityType;
//    const int dimDomain = DomainType::dimension;
//    
//    Eigen::MatrixXd perm(grid.size(0), dimDomain);
////    Eigen::VectorXd perm_x(grid.size(0));
////    Eigen::VectorXd perm_y(grid.size(0));
////    Eigen::VectorXd perm_z(grid.size(0));
//    JacobianRangeType permeability(0.0);
//
//    for (const auto& entity : DSC::viewRange(grid.leafView())) {
//      // get the barycenter
//      auto center = entity.geometry().center();
//      int index = grid.leafView().indexSet().index(entity);
//      for (const auto& dim : DSC::valueRange(dimDomain)) {
//        JacobianRangeType dir(0.0);
//        dir[0][dim] = 1.0;
//        diffusiveFlux(center, dir, permeability);
//        perm(index, dim) = permeability[0][dim];
//      }
////      perm_x(index) = permeabilityMatrix(0, 0);
////      perm_y(index) = permeabilityMatrix(1, 1);
////      perm_z(index) = permeabilityMatrix(2, 2);
//    }
//    
//    VTKWriter< typename CommonTraits::GridType::LeafGridView > vtkwriter(grid.leafView());
//    for (const auto& dim : DSC::valueRange(dimDomain))
//      vtkwriter.addCellData(perm.col(dim), "perm_" + DSC::toString(dim));
////    vtkwriter.addCellData(perm_y, "perm_y");
////    vtkwriter.addCellData(perm_z, "perm_z");
//    vtkwriter.write("data/Permeability");
//    
//    return;
//  }         // visualizePermeability

} // namespace SPE10
} // namespace Problem
} // namespace Multiscale {
} // namespace Dune {
