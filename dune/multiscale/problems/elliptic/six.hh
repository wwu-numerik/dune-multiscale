// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_SIX
#define DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_SIX

#include <dune/fem/function/common/function.hh>
#include <dune/multiscale/problems/base.hh>

#include <string>

#include "dune/multiscale/common/traits.hh"

namespace Dune {
namespace Multiscale {
namespace Problem {
/** \addtogroup problem_6 Problem::Six
 * @{ **/
//! ------------ Elliptic Problem 6 -------------------

namespace Six {

typedef CommonTraits::FunctionSpaceType::DomainType DomainType;
typedef CommonTraits::FunctionSpaceType::RangeType RangeType;
typedef CommonTraits::FunctionSpaceType::JacobianRangeType JacobianRangeType;
typedef CommonTraits::FunctionSpaceType::DomainFieldType TimeType;

struct ModelProblemData : public IModelProblemData {
  static const bool has_exact_solution = false;

  ModelProblemData();

  std::string getMacroGridFile() const;
  bool problemIsPeriodic() const;
  bool problemAllowsStochastics() const;
};

class Diffusion : public DiffusionBase {
public:
  void diffusiveFlux(const DomainType& x, const JacobianRangeType& gradient,
                     JacobianRangeType& flux) const;
  void jacobianDiffusiveFlux(const DomainType& x, const JacobianRangeType& /*position_gradient*/,
                             const JacobianRangeType& direction_gradient, JacobianRangeType& flux) const;
};

class ExactSolution : public Dune::Multiscale::CommonTraits::FunctionBaseType {
public:
  void evaluate(const DomainType& /*x*/, RangeType& /*y*/) const;
  void jacobian(const DomainType& /*x*/, JacobianRangeType& /*grad_u*/) const;
  void evaluate(const DomainType& x, const TimeType& /*timedummy*/, RangeType& y) const;
};

class DirichletData : public ZeroDirichletData {};
class NeumannData : public ZeroNeumannData {};
class LowerOrderTerm : public ZeroLowerOrder {};

MSCONSTANTFUNCTION(MassTerm, 0.0)
MSNULLFUNCTION(DirichletBoundaryCondition)
MSNULLFUNCTION(NeumannBoundaryCondition)
MSNULLFUNCTION(DefaultDummyFunction)
MSCONSTANTFUNCTION(FirstSource, 1.0)
MSNULLFUNCTION(SecondSource)

} // namespace Six {
}
} // namespace Multiscale {
} // namespace Dune {

#endif // ifndef DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_SIX
