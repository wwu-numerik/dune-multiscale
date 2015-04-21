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
  struct FunctionName##Mapper {                                                                                        \
    typedef std::function<ReturnType*()> FF;                                                                           \
    static std::map<std::string, FF> mk_map() {                                                                        \
      std::map<std::string, FF> funcs;                                                                                 \
      funcs.emplace("Synthetic", []() { return new DMP::Synthetic::FunctionName(); });                                 \
      funcs.emplace("Random", []() { return new DMP::Random::FunctionName(); });                                       \
      funcs.emplace("SPE10", []() { return new DMP::SPE10::FunctionName(); });                                         \
      funcs.emplace("Tarbert", []() { return new DMP::Tarbert::FunctionName(); });                                     \
      return funcs;                                                                                                    \
    }                                                                                                                  \
  };                                                                                                                   \
  static const auto FunctionName##_map = FunctionName##Mapper::mk_map();

FUNCTION_MAP(Problem::IModelProblemData, ModelProblemData)
FUNCTION_MAP(CommonTraits::FunctionBaseType, Source)
FUNCTION_MAP(CommonTraits::FunctionBaseType, ExactSolution)
FUNCTION_MAP(Problem::DiffusionBase, Diffusion)
FUNCTION_MAP(Problem::DirichletDataBase, DirichletData)
FUNCTION_MAP(Problem::NeumannDataBase, NeumannData)

template <class FunctionType>
FunctionType* make_f(const std::map<std::string, std::function<FunctionType*()>>& rets, std::string name) {
  const auto it = rets.find(name);
  if (it == rets.end())
    DUNE_THROW(Dune::InvalidStateException, "no data for Problem. (toggle PROBLEM_NINE_ONLY?)");
  const auto& el = it->second;
  return el();
}

Problem::ProblemContainer::ProblemContainer(MPIHelper::MPICommunicator global, MPIHelper::MPICommunicator local,
                                            DSC::Configuration config)
  : config_(config)
  , name_(config_.get("problem.name", "Synthetic"))
  , data_(make_f(ModelProblemData_map, name_))
  , source_(make_f(Source_map, name_))
  , exact_solution_(make_f(ExactSolution_map, name_))
  , diffusion_(make_f(Diffusion_map, name_))
  , dirichlet_(make_f(DirichletData_map, name_))
  , neumann_(make_f(NeumannData_map, name_)) {
  getMutableModelData().problem_init(*this, global, local);
}

const CommonTraits::FunctionBaseType& DMP::ProblemContainer::getSource() const { return *source_; }

const CommonTraits::FunctionBaseType& DMP::ProblemContainer::getExactSolution() const { return *exact_solution_; }

const Problem::IModelProblemData& DMP::ProblemContainer::getModelData() const { return *data_; }

Problem::IModelProblemData& DMP::ProblemContainer::getMutableModelData() { return *data_; }

const Problem::DiffusionBase& DMP::ProblemContainer::getDiffusion() const { return *diffusion_; }

Problem::DiffusionBase& DMP::ProblemContainer::getMutableDiffusion() {
  const auto& diffusion = getDiffusion();
  return const_cast<Problem::DiffusionBase&>(diffusion);
}

const Problem::DirichletDataBase& DMP::ProblemContainer::getDirichletData() const { return *dirichlet_; }

const Problem::NeumannDataBase& DMP::ProblemContainer::getNeumannData() const { return *neumann_; }

const std::string DMP::ProblemContainer::name() const { return name_; }

const Dune::Stuff::Common::Configuration& Problem::ProblemContainer::config() const { return config_; }

Dune::Stuff::Common::Configuration& Problem::ProblemContainer::config() { return config_; }
