// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_NINE
#define DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_NINE

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
/** \addtogroup problem_9 Problem::Nine
 * @{ **/
//! ------------ Elliptic Problem 9 -------------------

namespace Nine {

struct ModelProblemData : public IModelProblemData {
  virtual bool hasExactSolution() const { return true; }

  ModelProblemData();

  std::string getMacroGridFile() const;
  const BoundaryInfoType& boundaryInfo() const;
  const SubBoundaryInfoType& subBoundaryInfo() const;
  std::pair<CommonTraits::DomainType, CommonTraits::DomainType> gridCorners() const;

private:
  Dune::ParameterTree boundary_settings() const;
  std::unique_ptr<DSG::BoundaryInfos::NormalBased<typename View::Intersection>> boundaryInfo_;
  DSG::BoundaryInfos::AllDirichlet<typename SubView::Intersection> subBoundaryInfo_;
};

class Source : public Dune::Multiscale::CommonTraits::FunctionBaseType {
public:
  Source();

  PURE HOT void evaluate(const DomainType& x, RangeType& y) const;
  virtual size_t order() const { return 3; }
};

class Diffusion : public DiffusionBase {
public:
  Diffusion();

  //! currently used in gdt assembler
  virtual void evaluate(const DomainType& x, DiffusionBase::RangeType& y) const;
  PURE HOT void diffusiveFlux(const DomainType& x, const Problem::JacobianRangeType& direction,
                              Problem::JacobianRangeType& flux) const;

  //! HMM-only
  PURE void jacobianDiffusiveFlux(const DomainType& x, const Problem::JacobianRangeType& /*position_gradient*/,
                                  const Problem::JacobianRangeType& direction_gradient,
                                  Problem::JacobianRangeType& flux) const;
  virtual size_t order() const;
};

class ExactSolution : public Dune::Multiscale::CommonTraits::FunctionBaseType {
public:
  ExactSolution();

  PURE HOT void evaluate(const DomainType& x, RangeType& y) const;
  PURE HOT void jacobian(const DomainType& x, JacobianRangeType& grad_u) const;
  virtual size_t order() const { return 3; }
};

class DirichletData : public DirichletDataBase {
public:
  DirichletData() {}

  PURE void evaluate(const DomainType& x, RangeType& y) const;
  PURE void jacobian(const DomainType& x, JacobianRangeType& y) const;

private:
  ExactSolution solution_;
};

class NeumannData : public NeumannDataBase {
public:
  NeumannData() {}

  PURE void evaluate(const DomainType& x, RangeType& y) const;
};

MSNULLFUNCTION(DirichletBoundaryCondition)
MSNULLFUNCTION(NeumannBoundaryCondition)

} //! @} namespace Nine {
}
} // namespace Multiscale {
} // namespace Dune {

#endif // ifndef DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_NINE
