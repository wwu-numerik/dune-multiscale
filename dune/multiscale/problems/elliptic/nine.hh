// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_NINE
#define DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_NINE

#include <dune/fem/function/common/function.hh>
#include <dune/multiscale/problems/base.hh>

#include <memory>
#include <string>

#include "dune/multiscale/common/traits.hh"

namespace Dune {
namespace Multiscale {
namespace Problem {
/** \addtogroup problem_9 Problem::Nine
 * @{ **/
//! ------------ Elliptic Problem 9 -------------------

namespace Nine {

struct ModelProblemData : public IModelProblemData {
  virtual bool hasExactSolution() const { return true; }

  ModelProblemData();

  std::string getMacroGridFile() const;
  bool problemIsPeriodic() const;
  bool problemAllowsStochastics() const;
  std::unique_ptr<BoundaryInfoType> boundaryInfo() const;
  std::unique_ptr<SubBoundaryInfoType> subBoundaryInfo() const;
};

class FirstSource : public Dune::Multiscale::CommonTraits::FunctionBaseType {
public:
  FirstSource();

  void evaluate(const DomainType& x, RangeType& y) const;
  void evaluate(const DomainType& x, const TimeType& /*time*/, RangeType& y) const;
};

class Diffusion : public DiffusionBase {
public:
  Diffusion();

  void diffusiveFlux(const DomainType& x, const JacobianRangeType& direction, JacobianRangeType& flux) const;
  void jacobianDiffusiveFlux(const DomainType& x, const JacobianRangeType& /*position_gradient*/,
                             const JacobianRangeType& direction_gradient, JacobianRangeType& flux) const;
};

class ExactSolution : public Dune::Multiscale::CommonTraits::FunctionBaseType {
public:
  ExactSolution();

  void evaluate(const DomainType& x, RangeType& y) const;
  void jacobian(const DomainType& x, JacobianRangeType& grad_u) const;
};

class DirichletData : public DirichletDataBase {
public:
  DirichletData() {}

  void evaluate(const DomainType& x, RangeType& y) const;
  void evaluate(const DomainType& x, const TimeType& /*time*/, RangeType& y) const;
  void jacobian(const DomainType& x, JacobianRangeType& y) const;
  void jacobian(const DomainType& x, const TimeType& /*time*/, JacobianRangeType& y) const;
};

class NeumannData : public NeumannDataBase {
public:
  NeumannData() {}

  void evaluate(const DomainType& x, RangeType& y) const;
  void evaluate(const DomainType& x, const TimeType& /*time*/, RangeType& y) const;
};

class LowerOrderTerm : public ZeroLowerOrder {};

MSNULLFUNCTION(DirichletBoundaryCondition)
MSNULLFUNCTION(NeumannBoundaryCondition)
MSCONSTANTFUNCTION(MassTerm, 0.0)
MSNULLFUNCTION(DefaultDummyFunction)
MSNULLFUNCTION(SecondSource)

} //! @} namespace Nine {
}
} // namespace Multiscale {
} // namespace Dune {

#endif // ifndef DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_NINE
