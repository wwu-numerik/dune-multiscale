// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_MS_PROBLEMS_SELECTOR_HH
#define DUNE_MS_PROBLEMS_SELECTOR_HH

#include <dune/multiscale/msfem/msfem_traits.hh>
#include <dune/multiscale/problems/base.hh>
#include <memory>
#include <string>

#include "dune/multiscale/common/traits.hh"

namespace Dune {
template <class GridImp, class IntersectionImp>
class Intersection;
} // namespace Dune

namespace Dune {
namespace Multiscale {
namespace Problem {

typedef std::unique_ptr<const CommonTraits::FunctionBaseType> BasePtr;
const BasePtr& getFirstSource();
const BasePtr& getExactSolution();
const std::unique_ptr<const CommonTraits::DiffusionType>& getDiffusion();
const std::unique_ptr<const CommonTraits::DirichletBCType>& getDirichletBC();
const std::unique_ptr<const CommonTraits::NeumannBCType>& getNeumannBC();
const std::unique_ptr<const CommonTraits::ModelProblemDataType>& getModelData();
const std::unique_ptr<const CommonTraits::DirichletDataType>& getDirichletData();
const std::unique_ptr<const CommonTraits::NeumannDataType>& getNeumannData();

const std::string& name();

template <class IntersectionType>
typename std::enable_if<std::is_same<IntersectionType, CommonTraits::GridPartType::IntersectionType>::value, bool>::type
is_neumann(const IntersectionType& face) {
  return getModelData()->boundaryInfo().neumann(face);
}

template <class IntersectionType>
typename std::enable_if<!std::is_same<IntersectionType, CommonTraits::GridPartType::IntersectionType>::value,
                        bool>::type
is_neumann(const IntersectionType& face) {
  return getModelData()->subBoundaryInfo().neumann(face);
}

template <class IntersectionType>
typename std::enable_if<std::is_same<IntersectionType, CommonTraits::GridPartType::IntersectionType>::value, bool>::type
is_dirichlet(const IntersectionType& face) {
  return getModelData()->boundaryInfo().dirichlet(face);
}

template <class IntersectionType>
typename std::enable_if<!std::is_same<IntersectionType, CommonTraits::GridPartType::IntersectionType>::value,
                        bool>::type
is_dirichlet(const IntersectionType& face) {
  return getModelData()->subBoundaryInfo().dirichlet(face);
}

} //! @} namespace Problem
} // namespace Multiscale
} // namespace Dune

namespace DMP = Dune::Multiscale::Problem;

#endif // DUNE_MS_PROBLEMS_SELECTOR_HH
