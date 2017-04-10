// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_ER2007
#define DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_ER2007

#include <dune/multiscale/problems/base.hh>
#include <dune/stuff/grid/boundaryinfo.hh>
#include <memory>
#include <string>

#include "dune/multiscale/common/traits.hh"

#ifdef __GNUC__

#define PURE __attribute__((const))
#define HOT __attribute__((hot))
#define ALWAYS_INLINE __attribute__((always_inline)) inline

#else

#define PURE
#define HOT
#define ALWAYS_INLINE inline

#endif

namespace Dune {
namespace Multiscale {
namespace Problem {
/**
 *
------------ ER2007 Elliptic Problem -------------------
    out << "+============================================================+\n"
        << "|+==========================================================+|\n"
        << "||  Testcase ER2007: smooth data, nonhomogeneous dirichlet  ||\n"
        << "||  (see page 858 in Epshteyn, Riviere, 2007)               ||\n"
        << "|+----------------------------------------------------------+|\n"
        << "||  domain = [0, 1] x [0, 1]                                ||\n"
        << "||  diffusion = 1                                           ||\n"
        << "||  force     = 64 pi^2 (cos(8 pi x) + cos(8 pi y))         ||\n"
        << "||  dirichlet = cos(8 pi x) + cos(8 pi y)                   ||\n"
        << "||  exact solution = cos(8 pi x) + cos(8 pi y)              ||\n"
        << "|+==========================================================+|\n"
        << "+============================================================+" << std::endl;
**/

namespace ER2007 {

struct ModelProblemData : public IModelProblemData
{
  virtual bool hasExactSolution() const final override
  {
    return true;
  }

  ModelProblemData(MPIHelper::MPICommunicator /*global*/,
                   MPIHelper::MPICommunicator /*local*/,
                   Dune::XT::Common::Configuration /*config_in*/);

  const BoundaryInfoType& boundaryInfo() const final override;
  const SubBoundaryInfoType& subBoundaryInfo() const final override;
  std::pair<CommonTraits::DomainType, CommonTraits::DomainType> gridCorners() const final override;

private:
  std::unique_ptr<DSG::BoundaryInfos::AllDirichlet<typename View::Intersection>> boundaryInfo_;
  DSG::BoundaryInfos::AllDirichlet<typename SubView::Intersection> subBoundaryInfo_;
};

static const std::vector<std::string> exact_deriv{"-8.0*pi*sin(8.0*pi*x[0])", "-8.0*pi*sin(8.0*pi*x[1])"};

class Diffusion : public DiffusionBase
{
public:
  Diffusion(MPIHelper::MPICommunicator /*global*/,
            MPIHelper::MPICommunicator /*local*/,
            Dune::XT::Common::Configuration /*config_in*/)
  {
  }

  void diffusiveFlux(const DomainType& x,
                     const Problem::JacobianRangeType& direction,
                     Problem::JacobianRangeType& flux) const final override;
  void evaluate(const DomainType& x, RangeType& y) const final override;
};

MSEXPRESSIONFUNCTION(Source, "64*pi*pi*(cos(8.0*pi*x[0])+cos(8.0*pi*x[1]))", 3, exact_deriv)

using NeumannData = ZeroNeumannData;

MSEXPRESSIONFUNCTION(ExactSolution, "cos(8.0*pi*x[0])+cos(8.0*pi*x[1])", 3, exact_deriv)

class DirichletData : public DirichletDataBase
{
public:
  DirichletData(MPIHelper::MPICommunicator global,
                MPIHelper::MPICommunicator local,
                Dune::XT::Common::Configuration config_in)
    : solution_(global, local, config_in)
  {
  }

  PURE void evaluate(const DomainType& x, RangeType& y) const final override;
  PURE void jacobian(const DomainType& x, JacobianRangeType& y) const final override;

private:
  ExactSolution solution_;
};

} //! @} namespace ER2007 {
}
} // namespace Multiscale {
} // namespace Dune {

#endif // ifndef DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_ER2007
