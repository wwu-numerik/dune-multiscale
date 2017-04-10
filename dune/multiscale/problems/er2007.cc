#include <config.h>
#include <dune/common/exceptions.hh>
#include <dune/stuff/grid/boundaryinfo.hh>
#include <dune/xt/common/configuration.hh>
#include <dune/xt/common/validation.hh>
#include <math.h>
#include <sstream>

#include "dune/multiscale/problems/base.hh"
#include "er2007.hh"

namespace Dune {
namespace Multiscale {
namespace Problem {
namespace ER2007 {

ModelProblemData::ModelProblemData(MPIHelper::MPICommunicator global,
                                   MPIHelper::MPICommunicator local,
                                   Dune::XT::Common::Configuration config_in)
  : IModelProblemData(global, local, config_in)
  , boundaryInfo_(DSG::BoundaryInfos::AllDirichlet<typename View::Intersection>::create())
  , subBoundaryInfo_()
{
}

std::pair<CommonTraits::DomainType, CommonTraits::DomainType> ModelProblemData::gridCorners() const
{
  CommonTraits::DomainType lowerLeft(0.0);
  CommonTraits::DomainType upperRight(1.0);
  return {lowerLeft, upperRight};
}

const ModelProblemData::BoundaryInfoType& ModelProblemData::boundaryInfo() const
{
  return *boundaryInfo_;
}

const ModelProblemData::SubBoundaryInfoType& ModelProblemData::subBoundaryInfo() const
{
  return subBoundaryInfo_;
}

PURE void DirichletData::evaluate(const DomainType& x, RangeType& y) const
{
  solution_.evaluate(x, y);
} // evaluate

PURE void DirichletData::jacobian(const DomainType& x, JacobianRangeType& y) const
{
  solution_.jacobian(x, y);
}

void Diffusion::evaluate(const DomainType& x, Diffusion::RangeType& ret) const
{
  ret[0][0] = 1;
  ret[1][1] = 1;
}

void Diffusion::diffusiveFlux(const DomainType& x,
                              const Problem::JacobianRangeType& direction,
                              Problem::JacobianRangeType& flux) const
{
  flux[0][0] = direction[0][0];
  flux[0][1] = direction[0][1];
} // diffusiveFlux


} // namespace Nine
} // namespace Problem
} // namespace Multiscale {
} // namespace Dune {
