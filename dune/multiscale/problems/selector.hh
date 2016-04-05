// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_MS_PROBLEMS_SELECTOR_HH
#define DUNE_MS_PROBLEMS_SELECTOR_HH

#include <dune/multiscale/msfem/msfem_traits.hh>
#include <dune/multiscale/problems/base.hh>
#include <dune/xt/common/configuration.hh>
#include <memory>
#include <string>

#include "dune/multiscale/common/traits.hh"

namespace Dune {

template <class GridImp, class IntersectionImp>
class Intersection;

namespace Multiscale {
namespace Problem {

struct ProblemContainer {
  ProblemContainer(MPIHelper::MPICommunicator global, MPIHelper::MPICommunicator local, Dune::XT::Common::Configuration config_in);

  typedef std::unique_ptr<const CommonTraits::FunctionBaseType> BasePtr;
  const CommonTraits::FunctionBaseType& getSource() const;
  const CommonTraits::FunctionBaseType& getExactSolution() const;
  const DiffusionBase& getDiffusion() const;
  DiffusionBase& getMutableDiffusion();
  const IModelProblemData& getModelData() const;
  IModelProblemData& getMutableModelData();
  const DirichletDataBase& getDirichletData() const;
  const NeumannDataBase& getNeumannData() const;

  const std::string name() const;

  const Dune::XT::Common::Configuration& config() const;
  Dune::XT::Common::Configuration& config();

private:
  Dune::XT::Common::Configuration config_;
  const std::string name_;
  const std::unique_ptr<Problem::IModelProblemData> data_;
  const std::unique_ptr<const CommonTraits::FunctionBaseType> source_;
  const std::unique_ptr<const CommonTraits::FunctionBaseType> exact_solution_;
  const std::unique_ptr<Problem::DiffusionBase> diffusion_;
  const std::unique_ptr<const Problem::DirichletDataBase> dirichlet_;
  const std::unique_ptr<const Problem::NeumannDataBase> neumann_;
}; // struct ProblemContainer {

} //! @} namespace Problem
} // namespace Multiscale
} // namespace Dune

namespace DMP = Dune::Multiscale::Problem;

#endif // DUNE_MS_PROBLEMS_SELECTOR_HH
