#include <config.h>
#include <dune/common/exceptions.hh>
#include <dune/stuff/common/memory.hh>
#include <dune/stuff/common/parameter/configcontainer.hh>
#include <functional>
#include <map>

#include "dune/multiscale/problems/base.hh"
// for i in $(ls *hh) ; do echo \#include \"${i}\" ; done
#include "elliptic/nine.hh"
#include "elliptic/thirteen.hh"
#include "elliptic/spe10.hh"
#include "selector.hh"

namespace Dune {
template <class GridImp, class IntersectionImp>
class Intersection;
} // namespace Dune

#define PROBLEM_NAME Nine

using namespace Dune::Multiscale;

#define MAP_ITEM(ProblemName, ReturnType, FunctionName)                                                                \
  {                                                                                                                    \
    #ProblemName, []() { return static_cast<ReturnType>(DSC::make_unique<Problem::ProblemName::FunctionName>()); }     \
  }

#if PROBLEM_NINE_ONLY
# define FUNCTION_MAP(ReturnType, FunctionName)                                                                        \
    std::map<std::string, std::function<ReturnType()>>(                                                                \
         {MAP_ITEM(Nine, ReturnType, FunctionName), })
#else
# define FUNCTION_MAP(ReturnType, FunctionName)                                                                        \
    std::map<std::string, std::function<ReturnType()>>(                                                                \
        {MAP_ITEM(Nine, ReturnType, FunctionName), MAP_ITEM(SPE10, ReturnType, FunctionName),})
#endif

/* to add a new problem a line like this above
 * MAP_ITEM(NewProblemName, ReturnType, FunctionName), \
*/

// Toy doesn't comply with interface
// MAP_ITEM(Toy, ReturnType, FunctionName)

template <class FunctionType>
typename FunctionType::result_type find_and_call_item(const std::map<std::string, FunctionType>& rets) {
  auto it = rets.find(Dune::Multiscale::Problem::name());
  if (it == rets.end())
    DUNE_THROW(Dune::InvalidStateException, "no data for Problem. (toggle PROBLEM_NINE_ONLY?)");
  return it->second.operator()();
}

Problem::BasePtr Dune::Multiscale::Problem::getFirstSource() {
  static auto funcs = FUNCTION_MAP(BasePtr, FirstSource);
  return find_and_call_item(funcs);
}

Problem::BasePtr Dune::Multiscale::Problem::getSecondSource() {
  static auto funcs = FUNCTION_MAP(BasePtr, SecondSource);
  return find_and_call_item(funcs);
}

Problem::BasePtr Dune::Multiscale::Problem::getExactSolution() {
  static auto funcs = FUNCTION_MAP(BasePtr, ExactSolution);
  return find_and_call_item(funcs);
}

Problem::BasePtr Dune::Multiscale::Problem::getMassTerm() {
  static auto funcs = FUNCTION_MAP(BasePtr, MassTerm);
  return find_and_call_item(funcs);
}

Problem::BasePtr Dune::Multiscale::Problem::getDefaultDummyFunction() {
  static auto funcs = FUNCTION_MAP(BasePtr, DefaultDummyFunction);
  return find_and_call_item(funcs);
}

std::unique_ptr<const CommonTraits::ModelProblemDataType> Dune::Multiscale::Problem::getModelData() {
  static auto funcs = FUNCTION_MAP(std::unique_ptr<const CommonTraits::ModelProblemDataType>, ModelProblemData);
  return find_and_call_item(funcs);
}

std::unique_ptr<const CommonTraits::LowerOrderTermType> Dune::Multiscale::Problem::getLowerOrderTerm() {
  static auto funcs = FUNCTION_MAP(std::unique_ptr<const CommonTraits::LowerOrderTermType>, LowerOrderTerm);
  return find_and_call_item(funcs);
}

std::unique_ptr<const CommonTraits::DiffusionType> Dune::Multiscale::Problem::getDiffusion() {
  static auto funcs = FUNCTION_MAP(std::unique_ptr<const CommonTraits::DiffusionType>, Diffusion);
  return find_and_call_item(funcs);
}

std::unique_ptr<const CommonTraits::DirichletDataType> Dune::Multiscale::Problem::getDirichletData() {
  static auto funcs = FUNCTION_MAP(std::unique_ptr<const CommonTraits::DirichletDataType>, DirichletData);
  return find_and_call_item(funcs);
}

std::unique_ptr<const CommonTraits::NeumannDataType> Dune::Multiscale::Problem::getNeumannData() {
  static auto funcs = FUNCTION_MAP(std::unique_ptr<const CommonTraits::NeumannDataType>, NeumannData);
  return find_and_call_item(funcs);
}

std::unique_ptr<const CommonTraits::DirichletBCType> Dune::Multiscale::Problem::getDirichletBC() {
  static auto funcs = FUNCTION_MAP(std::unique_ptr<const CommonTraits::DirichletBCType>, DirichletBoundaryCondition);
  return find_and_call_item(funcs);
}

std::unique_ptr<const CommonTraits::NeumannBCType> Dune::Multiscale::Problem::getNeumannBC() {
  static auto funcs = FUNCTION_MAP(std::unique_ptr<const CommonTraits::NeumannBCType>, NeumannBoundaryCondition);
  return find_and_call_item(funcs);
}

std::string Dune::Multiscale::Problem::name() { return DSC_CONFIG_GET("problem.name", "Nine"); }
