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
  virtual bool hasExactSolution() const final override { return false; }

  ModelProblemData(MPIHelper::MPICommunicator /*global*/, MPIHelper::MPICommunicator /*local*/,
                   Dune::XT::Common::Configuration /*config_in*/);

  const BoundaryInfoType& boundaryInfo() const final override;
  const SubBoundaryInfoType& subBoundaryInfo() const final override;
  std::pair<CommonTraits::DomainType, CommonTraits::DomainType> gridCorners() const final override;

private:
  Dune::ParameterTree boundary_settings() const;
  std::unique_ptr<BoundaryInfoType> boundaryInfo_;
  std::unique_ptr<SubBoundaryInfoType> subBoundaryInfo_;
};

class Source : public Dune::Multiscale::CommonTraits::FunctionBaseType {
public:
  Source(MPIHelper::MPICommunicator /*global*/, MPIHelper::MPICommunicator /*local*/,
         Dune::XT::Common::Configuration /*config_in*/);

  void evaluate(const DomainType& x, RangeType& y) const final override;
  virtual size_t order() const final override { return 3; }
};

class Diffusion : public DiffusionBase {
public:
  Diffusion(MPIHelper::MPICommunicator /*global*/, MPIHelper::MPICommunicator /*local*/,
            Dune::XT::Common::Configuration /*config_in*/);
  ~Diffusion();

  //! currently used in gdt assembler
  virtual void evaluate(const DomainType& x, RangeType& y) const final override;

  void diffusiveFlux(const DomainType& x, const Problem::JacobianRangeType& direction,
                     Problem::JacobianRangeType& flux) const final override;

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
  DirichletData(MPIHelper::MPICommunicator /*global*/, MPIHelper::MPICommunicator /*local*/,
                Dune::XT::Common::Configuration /*config_in*/) {}

  void evaluate(const typename CommonTraits::DomainType& x, typename CommonTraits::RangeType& y) const final override;
};

class NeumannData : public NeumannDataBase {
private:
public:
public:
  NeumannData(MPIHelper::MPICommunicator /*global*/, MPIHelper::MPICommunicator /*local*/,
              Dune::XT::Common::Configuration /*config_in*/) {}

  void evaluate(const typename CommonTraits::DomainType& x, typename CommonTraits::RangeType& y) const final override;
};

MSNULLFUNCTION(ExactSolution)

} //! @} namespace Tarbert {
}
} // namespace Multiscale {
} // namespace Dune {

#endif // ifndef DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_TARBERT
