// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_ELEVEN
#define DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_ELEVEN

#include <dune/fem/function/common/function.hh>
#include <dune/multiscale/problems/base.hh>
#include <dune/multiscale/problems/constants.hh>
#include <string>

#include "dune/multiscale/common/traits.hh"

namespace Dune {
namespace Multiscale {
namespace Problem {
/** \addtogroup problem_11 Problem::Eleven
 * @{ **/
//! ------------ Elliptic Problem 11 -------------------

namespace Eleven {

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
private:
  typedef Dune::Multiscale::CommonTraits::FunctionSpaceType FunctionSpaceType;

public:
  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType RangeType;

  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

  static const int dimDomain = DomainType::dimension;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

  typedef typename FunctionSpaceType::DomainFieldType TimeType;

public:
  FirstSource() {}

  void evaluate(const DomainType& x, RangeType& y) const;
  void evaluate(const DomainType& x, const TimeType& /*time*/, RangeType& y) const;
};

MSNULLFUNCTION(SecondSource)

class Diffusion : public DiffusionBase {
  typedef Dune::Multiscale::CommonTraits::FunctionSpaceType FunctionSpaceType;
  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType RangeType;
  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

public:
  Diffusion();

  void diffusiveFlux(const DomainType& /*x*/, const JacobianRangeType& direction, JacobianRangeType& flux) const;
  void jacobianDiffusiveFlux(const DomainType& /*x*/, const JacobianRangeType& /*position_gradient*/,
                             const JacobianRangeType& direction_gradient, JacobianRangeType& flux) const;
};

MSCONSTANTFUNCTION(MassTerm, 0.0)
MSNULLFUNCTION(DefaultDummyFunction)

class LowerOrderTerm : public LowerOrderBase {
private:
  typedef Dune::Multiscale::CommonTraits::FunctionSpaceType FunctionSpaceType;

public:
  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType RangeType;
  typedef typename FunctionSpaceType::DomainFieldType TimeType;

public:
  LowerOrderTerm() {}

  void evaluate(const DomainType& x, RangeType& y) const;
  void evaluate(const DomainType& x, const TimeType& /*time*/, RangeType& y) const;
  void evaluate(const DomainType& x, const RangeType& position, const JacobianRangeType& direction_gradient,
                RangeType& y) const;
  void position_derivative(const DomainType& x, const RangeType& position, const JacobianRangeType& direction_gradient,
                           RangeType& y) const;
  void direction_derivative(const DomainType& x, const RangeType& position, const JacobianRangeType& direction_gradient,
                            JacobianRangeType& y) const;
};

class DirichletBoundaryCondition : public Dune::Multiscale::CommonTraits::FunctionBaseType {
private:
  typedef Dune::Multiscale::CommonTraits::FunctionSpaceType FunctionSpaceType;

public:
  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType RangeType;

  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

  static const int dimDomain = DomainType::dimension;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

  typedef typename FunctionSpaceType::DomainFieldType TimeType;

public:
  DirichletBoundaryCondition() {}

  void evaluate(const DomainType& x, RangeType& y) const;
  void evaluate(const DomainType& x, const TimeType& /*time*/, RangeType& y) const;
};

class NeumannBoundaryCondition : public Dune::Multiscale::CommonTraits::FunctionBaseType {
private:
  typedef Dune::Multiscale::CommonTraits::FunctionSpaceType FunctionSpaceType;

public:
  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType RangeType;

  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

  static const int dimDomain = DomainType::dimension;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

  typedef typename FunctionSpaceType::DomainFieldType TimeType;

public:
  NeumannBoundaryCondition() {}

  void evaluate(const DomainType& x, RangeType& y) const;
  void evaluate(const DomainType& x, const TimeType& /*time*/, RangeType& y) const;
};

class ExactSolution : public Dune::Multiscale::CommonTraits::FunctionBaseType {
public:
  typedef Dune::Multiscale::CommonTraits::FunctionSpaceType FunctionSpaceType;

public:
  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType RangeType;

  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

  typedef typename FunctionSpaceType::DomainFieldType TimeType;

public:
  ExactSolution();

  void evaluate(const DomainType& x, RangeType& y) const;
  void jacobian(const DomainType& x, typename FunctionSpaceType::JacobianRangeType& grad_u) const;
  void evaluate(const DomainType& x, const TimeType&, RangeType& y) const;
};

class DirichletData : public ZeroDirichletData {};
class NeumannData : public ZeroNeumannData {};

} //! @} namespace Eleven {
}
} // namespace Multiscale {
} // namespace Dune {

#endif // ifndef DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_ELEVEN
