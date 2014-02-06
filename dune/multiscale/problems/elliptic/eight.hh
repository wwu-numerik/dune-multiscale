// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_EIGHT
#define DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_EIGHT

#include <dune/fem/function/common/function.hh>
#include <dune/multiscale/common/traits.hh>
#include <dune/multiscale/problems/base.hh>

#include <string>

namespace Dune {
namespace Multiscale {
namespace Problem {
/** \addtogroup problem_8 Problem::Eight
 * @{ **/
//! ------------ Elliptic Problem 8 -------------------

namespace Eight {

struct ModelProblemData : public IModelProblemData {
  static const bool has_exact_solution = true;

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
  void evaluate(const DomainType& x, const TimeType& /*time*/, RangeType& y) const;
};

class Diffusion : public DiffusionBase {

  double additivePart(const DomainType& x, const int i, const int j) const; // additivePart

public:
  void diffusiveFlux(const DomainType& x, const JacobianRangeType& gradient, JacobianRangeType& flux) const;
  void jacobianDiffusiveFlux(const DomainType& x, const JacobianRangeType& position_gradient,
                             const JacobianRangeType& direction_gradient, JacobianRangeType& flux) const;
};

class ExactSolution : public Dune::Multiscale::CommonTraits::FunctionBaseType {
public:
  void evaluate(const DomainType& x, RangeType& y) const;
  void evaluate(const DomainType& x, const TimeType& /*timedummy*/, RangeType& y) const;
  void jacobian(const DomainType&, JacobianRangeType&) const;
};

class DirichletData : public ZeroDirichletData {};
class NeumannData : public ZeroNeumannData {};
class LowerOrderTerm : public ZeroLowerOrder {};

MSNULLFUNCTION(SecondSource)
MSCONSTANTFUNCTION(MassTerm, 0.0)
MSNULLFUNCTION(DefaultDummyFunction)
MSNULLFUNCTION(DirichletBoundaryCondition)
MSNULLFUNCTION(NeumannBoundaryCondition)

} //! @} namespace Eight {
}
} // namespace Multiscale {
} // namespace Dune {

#endif // ifndef DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_EIGHT
