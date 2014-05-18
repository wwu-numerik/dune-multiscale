// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_TARBERT
#define DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_TARBERT

#include <dune/multiscale/problems/base.hh>

#include <string>
#include <vector>

#include "dune/multiscale/common/traits.hh"

namespace Dune {
namespace Multiscale {
namespace Problem {
/** \addtogroup problem_tarbert Problem::Tarbert
 * @{ **/
//! ------------ Tarbert Problem -------------------

namespace Tarbert {

struct ModelProblemData : public IModelProblemData {
  virtual bool hasExactSolution() const { return false; }

  ModelProblemData();

  std::string getMacroGridFile() const;
  bool problemIsPeriodic() const;
  bool problemAllowsStochastics() const;
  const BoundaryInfoType& boundaryInfo() const;
  const SubBoundaryInfoType& subBoundaryInfo() const;
  std::pair<CommonTraits::DomainType, CommonTraits::DomainType> gridCorners() const;

private:
  Dune::ParameterTree boundary_settings() const;
  std::unique_ptr<BoundaryInfoType> boundaryInfo_;
  std::unique_ptr<SubBoundaryInfoType> subBoundaryInfo_;
};

class Source : public Dune::Multiscale::CommonTraits::FunctionBaseType {
public:
  Source();

  void evaluate(const DomainType& x, RangeType& y) const;
};

class Diffusion : public DiffusionBase {
public:
  Diffusion();
  ~Diffusion();

  //! currently used in gdt assembler
  virtual void evaluate(const DomainType& x, RangeType& y) const;

  void diffusiveFlux(const DomainType& x, const Problem::JacobianRangeType& direction, Problem::JacobianRangeType& flux) const;
  void jacobianDiffusiveFlux(const DomainType& x, const Problem::JacobianRangeType& /*position_gradient*/,
                             const Problem::JacobianRangeType& direction_gradient, Problem::JacobianRangeType& flux) const;
  
//  void visualizePermeability(const CommonTraits::GridType& grid) const;
private:
  void readPermeability();

  std::vector<double> deltas_;
  double* permeability_; //! \todo automatic memory
  mutable DomainType permIntervalls_;
  mutable Dune::FieldMatrix<double, DomainType::dimension, DomainType::dimension> permMatrix_;
};

class DirichletData : public DirichletDataBase {
private:
public:
public:
  DirichletData() {}

  void evaluate(const typename FunctionSpaceType::DomainType& x, typename FunctionSpaceType::RangeType& y) const;
};

class NeumannData : public NeumannDataBase {
private:
public:
public:
  NeumannData() {}

  void evaluate(const typename FunctionSpaceType::DomainType& x, typename FunctionSpaceType::RangeType& y) const;
};


MSNULLFUNCTION(DirichletBoundaryCondition)
MSNULLFUNCTION(NeumannBoundaryCondition)
MSNULLFUNCTION(ExactSolution)

} //! @} namespace Tarbert {
}
} // namespace Multiscale {
} // namespace Dune {

#endif // ifndef DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_TARBERT
