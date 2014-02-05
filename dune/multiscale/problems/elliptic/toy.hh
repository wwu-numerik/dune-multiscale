// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_TOY
#define DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_TOY

#include <dune/fem/function/common/function.hh>
#include <dune/multiscale/problems/base.hh>
#include <dune/multiscale/problems/constants.hh>
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
  static const bool has_exact_solution = true;

  ModelProblemData();

  std::string getMacroGridFile() const;
  bool problemIsPeriodic() const;
  bool problemAllowsStochastics() const;
};

class FirstSource : public Dune::Multiscale::CommonTraits::FunctionBaseType {
private:
  typedef Dune::Multiscale::CommonTraits::FunctionSpaceType FunctionSpaceType;

public:
  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType RangeType;
  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

  static const int dimDomain = DomainType::dimension;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

  typedef DomainFieldType TimeType;

public:
  FirstSource() {}

  void evaluate(const DomainType& x, RangeType& y) const {
    double a_0_x_0 = 1.0 + pow(x[0], 2.0);
    double a_1_x_1 = 1.0 + pow(x[0], 2.0);

    double grad_a_0_x_0 = 2.0 * x[0];
    double grad_a_1_x_1 = 0.0;

    JacobianRangeType grad_u(0.0);

    grad_u[0][0] = (1.0 - x[0]) * (1.0 - x[1]) * x[1] - x[0] * (1.0 - x[1]) * x[1];
    grad_u[0][1] = x[0] * (1.0 - x[0]) * (1.0 - x[1]) - x[0] * (1.0 - x[0]) * x[1];

    const RangeType d_xx_u = (-2.0) * (1.0 - x[1]) * x[1];
    const RangeType d_yy_u = (-2.0) * (1.0 - x[0]) * x[0];

    y = 0.0;
    y -= grad_a_0_x_0 * grad_u[0][0];
    y -= a_0_x_0 * d_xx_u;

    y -= grad_a_1_x_1 * grad_u[0][1];
    y -= a_1_x_1 * d_yy_u;
  } // evaluate

  void evaluate(const DomainType& x, const TimeType& /*time*/, RangeType& y) const { evaluate(x, y); }
};

MSNULLFUNCTION(SecondSource)

class Diffusion : public DiffusionBase {
public:
  typedef Dune::Multiscale::CommonTraits::FunctionSpaceType FunctionSpaceType;

public:
  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType RangeType;
  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

  typedef DomainFieldType TimeType;

public:
  Diffusion() {}

  void diffusiveFlux(const DomainType& x, const JacobianRangeType& gradient, JacobianRangeType& flux) const;
};

class LowerOrderTerm : public ZeroLowerOrder {};

MSCONSTANTFUNCTION(MassTerm, 0.0)
MSNULLFUNCTION(DirichletBoundaryCondition)
MSNULLFUNCTION(NeumannBoundaryCondition)
MSNULLFUNCTION(DefaultDummyFunction)

class ExactSolution : public Dune::Multiscale::CommonTraits::FunctionBaseType {
public:
  typedef Dune::Multiscale::CommonTraits::FunctionSpaceType FunctionSpaceType;

public:
  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType RangeType;
  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

  typedef DomainFieldType TimeType;

public:
  ExactSolution() {}

  void evaluate(const DomainType& x, RangeType& y) const;
  void evaluate(const DomainType& x, const TimeType& /*timedummy*/, RangeType& y) const;
  void jacobian(const DomainType& x, JacobianRangeType& grad_u) const;
};

} //! @} namespace Toy {
}
} // namespace Multiscale {
} // namespace Dune {

#endif // ifndef DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_TOY
