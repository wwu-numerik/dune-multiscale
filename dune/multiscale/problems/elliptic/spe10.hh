// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_SPE10
#define DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_SPE10

#include <dune/fem/function/common/function.hh>
#include <dune/multiscale/problems/base.hh>

#include <string>
#include <vector>

#include "dune/multiscale/common/traits.hh"

namespace Dune {
namespace Multiscale {
namespace Problem {
/** \addtogroup problem_spe10 Problem::SPE10
 * @{ **/
//! ------------ SPE10 Problem -------------------

namespace SPE10 {

typedef CommonTraits::FunctionSpaceType::DomainType DomainType;
typedef CommonTraits::FunctionSpaceType::RangeType RangeType;
typedef CommonTraits::FunctionSpaceType::JacobianRangeType JacobianRangeType;
typedef CommonTraits::FunctionSpaceType::DomainFieldType TimeType;

struct ModelProblemData : public IModelProblemData {
  virtual bool hasExactSolution() const { return false; }

  ModelProblemData();

  std::string getMacroGridFile() const;
  bool problemIsPeriodic() const;
  bool problemAllowsStochastics() const;
  std::unique_ptr<BoundaryInfoType> boundaryInfo() const;
  std::unique_ptr<SubBoundaryInfoType> subBoundaryInfo() const;

private:
  Dune::ParameterTree boundary_settings() const;
};

class FirstSource : public Dune::Multiscale::CommonTraits::FunctionBaseType {
public:
  FirstSource();

  void evaluate(const DomainType& x, RangeType& y) const;
  void evaluate(const DomainType& x, const TimeType& /*time*/, RangeType& y) const;
};

MSNULLFUNCTION(SecondSource)

class Diffusion : public DiffusionBase {
public:
  Diffusion();
  ~Diffusion();

  void diffusiveFlux(const DomainType& x, const JacobianRangeType& direction, JacobianRangeType& flux) const;
  void jacobianDiffusiveFlux(const DomainType& x, const JacobianRangeType& /*position_gradient*/,
                             const JacobianRangeType& direction_gradient, JacobianRangeType& flux) const;

private:
  void readPermeability();

  std::vector<double> deltas_;
  double* permeability_;//! TODO automatic memory
};

class LowerOrderTerm : public ZeroLowerOrder {};

MSCONSTANTFUNCTION(MassTerm, 0.0)
MSNULLFUNCTION(DirichletBoundaryCondition)
MSNULLFUNCTION(NeumannBoundaryCondition)
MSNULLFUNCTION(DefaultDummyFunction)
MSNULLFUNCTION(ExactSolution)

class DirichletData : public DirichletDataBase {
private:


public:








public:
  DirichletData() {}

  void evaluate(const typename FunctionSpaceType::DomainType &x, typename FunctionSpaceType::RangeType &y) const;
  void evaluate(const typename FunctionSpaceType::DomainType &x, const TimeType &/*time*/,
                typename FunctionSpaceType::RangeType &y) const;
};

class NeumannData : public NeumannDataBase {
private:


public:








public:
  NeumannData() {}

  void evaluate(const typename FunctionSpaceType::DomainType &x, typename FunctionSpaceType::RangeType &y) const;
  void evaluate(const typename FunctionSpaceType::DomainType &x, const TimeType &/*time*/,
                typename FunctionSpaceType::RangeType &y) const;
};

} //! @} namespace SPE10 {
}
} // namespace Multiscale {
} // namespace Dune {

#endif // ifndef DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_NINE
