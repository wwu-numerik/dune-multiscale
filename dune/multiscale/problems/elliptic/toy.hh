// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_TOY
#define DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_TOY

#include <dune/fem/function/common/function.hh>
#include <dune/multiscale/problems/base.hh>

#include <string>

#include "dune/multiscale/common/traits.hh"

namespace Dune {
namespace Multiscale {
namespace Problem {
/** \addtogroup problem_0 Problem::Toy
 * @{ **/
// --------- TOY PROBLEM FOR SIMPLE TESTS --------------
//! ------------ Elliptic Problem 0 --------------------
// an epsilon-indepent linear elliptic toy model problem

// we solve for
//  u(x) = x_1 ( 1 - x_1 ) (1 - x_2 ) x_2
// with
//  a_12(x) = a_21(x) = 0 and a_11(x) = a_22(x) = 1 + (x_1)²
// f is defined so that everything fits together

namespace Toy {

struct ModelProblemData : public IModelProblemData {
  ModelProblemData();

  std::string getMacroGridFile() const;
  bool problemIsPeriodic() const;
  bool problemAllowsStochastics() const;
  bool hasExactSolution() const { return true; }
};

class FirstSource : public Dune::Multiscale::CommonTraits::FunctionBaseType {
public:
  void evaluate(const DomainType& x, RangeType& y) const; // evaluate
  void evaluate(const DomainType& x, const TimeType& /*time*/, RangeType& y) const { evaluate(x, y); }
};

class Diffusion : public DiffusionBase {
  void diffusiveFlux(const DomainType& x, const JacobianRangeType& gradient, JacobianRangeType& flux) const;
};

class ExactSolution : public Dune::Multiscale::CommonTraits::FunctionBaseType {
public:
  void evaluate(const DomainType& x, RangeType& y) const;
  void evaluate(const DomainType& x, const TimeType& /*timedummy*/, RangeType& y) const;
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

class LowerOrderTerm : public ZeroLowerOrder {};
class NeumannData : public ZeroNeumannData {};

MSCONSTANTFUNCTION(MassTerm, 0.0)
MSNULLFUNCTION(DirichletBoundaryCondition)
MSNULLFUNCTION(NeumannBoundaryCondition)
MSNULLFUNCTION(DefaultDummyFunction)
MSNULLFUNCTION(SecondSource)

} //! @} namespace Toy {
}
} // namespace Multiscale {
} // namespace Dune {

#endif // ifndef DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_TOY
