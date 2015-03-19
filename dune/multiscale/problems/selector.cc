#include <config.h>
#include <dune/common/exceptions.hh>
#include <dune/stuff/common/memory.hh>
#include <dune/stuff/common/parallel/threadstorage.hh>

#include <functional>
#include <unordered_map>
#include <map>

#include "dune/multiscale/problems/base.hh"
// for i in $(ls *hh) ; do echo \#include \"${i}\" ; done
#include "synthetic.hh"
#include "thirteen.hh"
#include "spe10.hh"
#include "tarbert.hh"
#include "random.hh"
#include "selector.hh"

namespace Dune {
template <class GridImp, class IntersectionImp>
class Intersection;
} // namespace Dune

using namespace Dune::Multiscale;

/* to add a new problem add another emplace line below
 * funcs.emplace("NAME", std::unique_ptr<const ReturnType>(new DMP::NAME::FunctionName())); \
*/
#define FUNCTION_MAP(ReturnType, FunctionName)                                                                         \
  struct foo {                                                                                                         \
    static std::map<std::string, const std::unique_ptr<const ReturnType>> mk_map() {                                   \
      std::map<std::string, const std::unique_ptr<const ReturnType>> funcs;                                            \
      funcs.emplace("Synthetic", std::unique_ptr<const ReturnType>(new DMP::Synthetic::FunctionName()));               \
      funcs.emplace("Random", std::unique_ptr<const ReturnType>(new DMP::Random::FunctionName()));                     \
      funcs.emplace("SPE10", std::unique_ptr<const ReturnType>(new DMP::SPE10::FunctionName()));                       \
      funcs.emplace("Tarbert", std::unique_ptr<const ReturnType>(new DMP::Tarbert::FunctionName()));                   \
      return funcs;                                                                                                    \
    }                                                                                                                  \
  };                                                                                                                   \
  static const auto funcs = foo::mk_map()

template <class FunctionType>
const FunctionType& find_and_call_item(const std::map<std::string, const std::unique_ptr<const FunctionType>>& rets) {
  const auto it = rets.find(Dune::Multiscale::Problem::name());
  if (it == rets.end())
    DUNE_THROW(Dune::InvalidStateException, "no data for Problem. (toggle PROBLEM_NINE_ONLY?)");
  const auto& el = it->second;
  return *el;
}

const CommonTraits::FunctionBaseType& Dune::Multiscale::Problem::getSource() {
  FUNCTION_MAP(CommonTraits::FunctionBaseType, Source);
  return find_and_call_item(funcs);
}

const CommonTraits::FunctionBaseType& Dune::Multiscale::Problem::getExactSolution() {
  FUNCTION_MAP(CommonTraits::FunctionBaseType, ExactSolution);
  return find_and_call_item(funcs);
}

const Problem::IModelProblemData& Dune::Multiscale::Problem::getModelData() {
  FUNCTION_MAP(Problem::IModelProblemData, ModelProblemData);
  return find_and_call_item(funcs);
}

Problem::IModelProblemData& Dune::Multiscale::Problem::getMutableModelData() {
  const auto& data = getModelData();
  return const_cast<Problem::IModelProblemData&>(data);
}

const Problem::DiffusionBase& Dune::Multiscale::Problem::getDiffusion() {
  FUNCTION_MAP(Problem::DiffusionBase, Diffusion);
  return find_and_call_item(funcs);
}

Problem::DiffusionBase& Dune::Multiscale::Problem::getMutableDiffusion() {
  const auto& diffusion = getDiffusion();
  return const_cast<Problem::DiffusionBase&>(diffusion);
}

const Problem::DirichletDataBase& Dune::Multiscale::Problem::getDirichletData() {
  FUNCTION_MAP(Problem::DirichletDataBase, DirichletData);
  return find_and_call_item(funcs);
}

const Problem::NeumannDataBase& Dune::Multiscale::Problem::getNeumannData() {
  FUNCTION_MAP(Problem::NeumannDataBase, NeumannData);
  return find_and_call_item(funcs);
}

const std::string& Dune::Multiscale::Problem::name() {
  static std::string myName = DSC_CONFIG_GET("problem.name", "Nine");
  return myName;
}
