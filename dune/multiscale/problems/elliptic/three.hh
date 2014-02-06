// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_THREE
#define DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_THREE

#include <dune/fem/function/common/function.hh>
#include <dune/multiscale/problems/base.hh>

#include <string>

#include "dune/multiscale/common/traits.hh"

namespace Dune {
namespace Multiscale {
namespace Problem {
/** \addtogroup problem_3 Problem::Three
 * @{ **/
//! ------------ Elliptic Problem 3 -------------------
// nonlinear elliptic model problem - heterogeneous setting
// no exact solution available!

namespace Three {

struct ModelProblemData : public IModelProblemData {
  static const bool has_exact_solution = false;

  ModelProblemData();

  std::string getMacroGridFile() const;
  bool problemIsPeriodic() const;
  bool problemAllowsStochastics() const;
  bool symmetricDiffusion() const { return false; }
  bool linear() const { return false; }
};

class FirstSource : public Dune::Multiscale::CommonTraits::FunctionBaseType {
public:
  void evaluate(const DomainType& x, RangeType& y) const;
  void evaluate(const DomainType& x, const TimeType& /*time*/, RangeType& y) const { evaluate(x, y); }
};

class Diffusion : public DiffusionBase {
public:
  void diffusiveFlux(const DomainType& x, const JacobianRangeType& gradient, JacobianRangeType& flux) const;
  void jacobianDiffusiveFlux(const DomainType& x, const JacobianRangeType& position_gradient,
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
MSNULLFUNCTION(SecondSource)

} //! @} namespace Three {
}
} // namespace Multiscale {
} // namespace Dune {

#endif // ifndef DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_THREE
