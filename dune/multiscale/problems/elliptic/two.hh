// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_TWO
#define DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_TWO

#include <dune/fem/function/common/function.hh>
#include <dune/multiscale/common/traits.hh>
#include <dune/multiscale/problems/base.hh>

#include <string>

namespace Dune {
namespace Multiscale {
namespace Problem {

/** \addtogroup problem_2 Problem::Two
 * @{ **/
//! ------------ Elliptic Problem 2 -------------------

namespace Two {

struct ModelProblemData : public IModelProblemData {
  static const bool has_exact_solution = false;

  ModelProblemData();

  std::string getMacroGridFile() const;
  bool problemIsPeriodic() const;
  bool problemAllowsStochastics() const;
  bool symmetricDiffusion() const { return false; }
};

class Diffusion : public DiffusionBase {
public:
  void diffusiveFlux(const DomainType& x, const JacobianRangeType& gradient, JacobianRangeType& flux) const;
  void jacobianDiffusiveFlux(const DomainType& x, const JacobianRangeType& /*position_gradient*/,
                             const JacobianRangeType& direction_gradient, JacobianRangeType& flux) const;
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
MSNULLFUNCTION(ExactSolution)

} //! @} namespace Two {
}
} // namespace Multiscale {
} // namespace Dune {
//! @}

#endif // ifndef DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_TWO
