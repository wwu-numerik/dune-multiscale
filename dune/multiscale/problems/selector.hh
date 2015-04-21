// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_MS_PROBLEMS_SELECTOR_HH
#define DUNE_MS_PROBLEMS_SELECTOR_HH

#include <dune/multiscale/msfem/msfem_traits.hh>
#include <dune/multiscale/problems/base.hh>
#include <dune/stuff/common/configuration.hh>
#include <memory>
#include <string>

#include "dune/multiscale/common/traits.hh"

namespace Dune {

template <class GridImp, class IntersectionImp>
class Intersection;

namespace Multiscale {
namespace Problem {

struct ProblemContainer {
  ProblemContainer(MPIHelper::MPICommunicator global, MPIHelper::MPICommunicator local, DSC::Configuration config);

  typedef std::unique_ptr<const CommonTraits::FunctionBaseType> BasePtr;
  const CommonTraits::FunctionBaseType& getSource();
  const CommonTraits::FunctionBaseType& getExactSolution();
  const DiffusionBase& getDiffusion();
  DiffusionBase& getMutableDiffusion();
  const IModelProblemData& getModelData();
  IModelProblemData& getMutableModelData();
  const DirichletDataBase& getDirichletData();
  const NeumannDataBase& getNeumannData();

  const std::string& name();

  template <class IntersectionType>
  typename std::enable_if<std::is_same<IntersectionType, CommonTraits::GridViewType::Intersection>::value, bool>::type
  is_neumann(const IntersectionType& face) {
    return getModelData().boundaryInfo().neumann(face);
  }

  template <class IntersectionType>
  typename std::enable_if<!std::is_same<IntersectionType, CommonTraits::GridViewType::Intersection>::value, bool>::type
  is_neumann(const IntersectionType& face) {
    return getModelData().subBoundaryInfo().neumann(face);
  }

  template <class IntersectionType>
  typename std::enable_if<std::is_same<IntersectionType, CommonTraits::GridViewType::Intersection>::value, bool>::type
  is_dirichlet(const IntersectionType& face) {
    return getModelData().boundaryInfo().dirichlet(face);
  }

  template <class IntersectionType>
  typename std::enable_if<!std::is_same<IntersectionType, CommonTraits::GridViewType::Intersection>::value, bool>::type
  is_dirichlet(const IntersectionType& face) {
    return getModelData().subBoundaryInfo().dirichlet(face);
  }
private:
  DSC::Configuration config_;
  const std::string name_;
  std::unique_ptr<Problem::IModelProblemData> data_;
  std::unique_ptr<CommonTraits::FunctionBaseType> source_;
  std::unique_ptr<CommonTraits::FunctionBaseType> exact_solution_;
  std::unique_ptr<Problem::DiffusionBase> diffusion_;
  std::unique_ptr<Problem::DirichletDataBase> dirichlet_;
  std::unique_ptr<Problem::NeumannDataBase> neumann_;
}; // struct ProblemContainer {

} //! @} namespace Problem
} // namespace Multiscale
} // namespace Dune

namespace DMP = Dune::Multiscale::Problem;

#endif // DUNE_MS_PROBLEMS_SELECTOR_HH
