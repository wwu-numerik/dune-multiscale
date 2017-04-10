// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_THIRTEEN
#define DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_THIRTEEN

#include <dune/multiscale/problems/base.hh>

#include <string>

#include "dune/multiscale/common/traits.hh"

namespace Dune {
namespace Multiscale {
namespace Problem {
/** \addtogroup problem_13 Problem::Thirteen
 * @{ **/
//! ------------ Elliptic Problem 13 -------------------

namespace Thirteen {

struct ModelProblemData : public IModelProblemData
{
  ModelProblemData();
};

class Source : public Dune::Multiscale::CommonTraits::FunctionBaseType
{
public:
  Source(MPIHelper::MPICommunicator /*global*/,
         MPIHelper::MPICommunicator /*local*/,
         Dune::XT::Common::Configuration /*config_in*/);

  void evaluate(const DomainType& x, RangeType& y) const final override;
};

class Diffusion : public DiffusionBase
{
public:
  Diffusion(MPIHelper::MPICommunicator /*global*/,
            MPIHelper::MPICommunicator /*local*/,
            Dune::XT::Common::Configuration /*config_in*/);

  void diffusiveFlux(const DomainType& x,
                     const Problem::JacobianRangeType& direction,
                     Problem::JacobianRangeType& flux) const final override;
};

class NeumannBoundaryCondition : public Dune::Multiscale::CommonTraits::FunctionBaseType
{
public:
  void evaluate(const DomainType& x, RangeType& y) const final override;
};

class DirichletData : public ZeroDirichletData
{
};
class NeumannData : public ZeroNeumannData
{
};

MSNULLFUNCTION(ExactSolution)

} //! @} namespace Thirteen {
}
} // namespace Multiscale {
} // namespace Dune {

#endif // ifndef DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_THIRTEEN
