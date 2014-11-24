#include <config.h>
#include <dune/common/exceptions.hh>
#include <dune/stuff/common/memory.hh>
#include <dune/stuff/common/parallel/threadstorage.hh>

#include <functional>
#include <unordered_map>
#include <map>

#include "dune/multiscale/problems/base.hh"
// for i in $(ls *hh) ; do echo \#include \"${i}\" ; done
#include "elliptic/nine.hh"
#include "elliptic/thirteen.hh"
#include "elliptic/spe10.hh"
#include "elliptic/tarbert.hh"
#include "selector.hh"

namespace Dune {
template <class GridImp, class IntersectionImp>
class Intersection;
} // namespace Dune

using namespace Dune::Multiscale;

template <class ReturnType>
struct AutoInitBase {

  virtual ~AutoInitBase() {}

  AutoInitBase()
    : ptr(nullptr) {}

  virtual ReturnType make() const = 0;
  const ReturnType& call() {
    if (ptr == nullptr)
      ptr = make();
    return ptr;
  }

  ReturnType ptr;
};

template <class ReturnType, class Function>
struct AutoInit : public AutoInitBase<ReturnType> {

  AutoInit()
    : AutoInitBase<ReturnType>()
    , maker([]() { return static_cast<ReturnType>(DSC::make_unique<Function>()); }) {}

  virtual ReturnType make() const { return maker.operator()(); }
  const std::function<ReturnType()> maker;
};

#define MAP_ITEM(ProblemName, ReturnType, FunctionName)                                                                \
  {                                                                                                                    \
    #ProblemName, std::shared_ptr < AutoInitBase < ReturnType >> (new AutoInit < ReturnType,                           \
                                                                  Problem::ProblemName::FunctionName > ())             \
  }

#define FUNCTION_MAP(ReturnType, FunctionName)                                                                         \
  DS::PerThreadValue< std::map<std::string, std::shared_ptr<AutoInitBase<ReturnType>>>> funcs({MAP_ITEM(Nine, ReturnType, FunctionName),          \
                                                                    MAP_ITEM(SPE10, ReturnType, FunctionName),         \
                                                                    MAP_ITEM(Tarbert, ReturnType, FunctionName)})

/* to add a new problem a line like this above
 * MAP_ITEM(NewProblemName, ReturnType, FunctionName), \
*/

template <class FunctionType>
const FunctionType& find_and_call_item(const std::map<std::string, std::shared_ptr<AutoInitBase<FunctionType>>>& rets) {
  auto it = rets.find(Dune::Multiscale::Problem::name());
  if (it == rets.end())
    DUNE_THROW(Dune::InvalidStateException, "no data for Problem. (toggle PROBLEM_NINE_ONLY?)");
  const auto& el = it->second->call();
  return el;
}

const Problem::BasePtr& Dune::Multiscale::Problem::getSource() {
  static const FUNCTION_MAP(BasePtr, Source);
  return find_and_call_item(*funcs);
}

const Problem::BasePtr& Dune::Multiscale::Problem::getExactSolution() {
  static const FUNCTION_MAP(BasePtr, ExactSolution);
  return find_and_call_item(*funcs);
}

const std::unique_ptr<const Problem::IModelProblemData>& Dune::Multiscale::Problem::getModelData() {
  static const FUNCTION_MAP(std::unique_ptr<const Problem::IModelProblemData>, ModelProblemData);
  return find_and_call_item(*funcs);
}

const std::unique_ptr<const Problem::DiffusionBase>& Dune::Multiscale::Problem::getDiffusion() {
  static const FUNCTION_MAP(std::unique_ptr<const Problem::DiffusionBase>, Diffusion);
  return find_and_call_item(*funcs);
}

const std::unique_ptr<const Problem::DirichletDataBase>& Dune::Multiscale::Problem::getDirichletData() {
  static const FUNCTION_MAP(std::unique_ptr<const Problem::DirichletDataBase>, DirichletData);
  return find_and_call_item(*funcs);
}

const std::unique_ptr<const Problem::NeumannDataBase>& Dune::Multiscale::Problem::getNeumannData() {
  static const FUNCTION_MAP(std::unique_ptr<const Problem::NeumannDataBase>, NeumannData);
  return find_and_call_item(*funcs);
}

const std::string& Dune::Multiscale::Problem::name() {
  static std::string myName = DSC_CONFIG_GET("problem.name", "Nine");
  return myName;
}
