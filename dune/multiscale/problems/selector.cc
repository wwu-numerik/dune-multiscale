#include <config.h>
#include <dune/common/exceptions.hh>
#include <dune/xt/common/memory.hh>
#include <dune/xt/common/misc.hh>
#include <dune/xt/common/parallel/threadstorage.hh>

#include <functional>
#include <map>
#include <unordered_map>

#include "dune/multiscale/problems/base.hh"
// for i in $(ls *hh) ; do echo \#include \"${i}\" ; done
#include "er2007.hh"
#include "random.hh"
#include "selector.hh"
#include "spe10.hh"
#include "synthetic.hh"
#include "tarbert.hh"
#include "thirteen.hh"

namespace Dune {
template <class GridImp, class IntersectionImp>
class Intersection;
} // namespace Dune

using namespace Dune::Multiscale;

/* to add a new problem add another emplace line below
 * Dune::XT::Common::map_emplace(funcs, "NAME", std::unique_ptr<const
 * ReturnType>(new DMP::NAME::FunctionName())); \
*/
#define FUNCTION_MAP(ReturnType, FunctionName)                                                                         \
  struct FunctionName##Mapper                                                                                          \
  {                                                                                                                    \
    typedef std::function<ReturnType*(                                                                                 \
        Dune::MPIHelper::MPICommunicator, Dune::MPIHelper::MPICommunicator, Dune::XT::Common::Configuration)>          \
        FF;                                                                                                            \
    static std::map<std::string, FF> mk_map()                                                                          \
    {                                                                                                                  \
      std::map<std::string, FF> funcs;                                                                                 \
      Dune::XT::Common::map_emplace(funcs,                                                                             \
                                    "Synthetic",                                                                       \
                                    [](Dune::MPIHelper::MPICommunicator global,                                        \
                                       Dune::MPIHelper::MPICommunicator local,                                         \
                                       Dune::XT::Common::Configuration config_in) {                                    \
                                      return new DMP::Synthetic::FunctionName(global, local, config_in);               \
                                    });                                                                                \
      Dune::XT::Common::map_emplace(funcs,                                                                             \
                                    "Random",                                                                          \
                                    [](Dune::MPIHelper::MPICommunicator global,                                        \
                                       Dune::MPIHelper::MPICommunicator local,                                         \
                                       Dune::XT::Common::Configuration config_in) {                                    \
                                      return new DMP::Random::FunctionName(global, local, config_in);                  \
                                    });                                                                                \
      Dune::XT::Common::map_emplace(funcs,                                                                             \
                                    "SPE10",                                                                           \
                                    [](Dune::MPIHelper::MPICommunicator global,                                        \
                                       Dune::MPIHelper::MPICommunicator local,                                         \
                                       Dune::XT::Common::Configuration config_in) {                                    \
                                      return new DMP::SPE10::FunctionName(global, local, config_in);                   \
                                    });                                                                                \
      Dune::XT::Common::map_emplace(funcs,                                                                             \
                                    "Tarbert",                                                                         \
                                    [](Dune::MPIHelper::MPICommunicator global,                                        \
                                       Dune::MPIHelper::MPICommunicator local,                                         \
                                       Dune::XT::Common::Configuration config_in) {                                    \
                                      return new DMP::Tarbert::FunctionName(global, local, config_in);                 \
                                    });                                                                                \
      Dune::XT::Common::map_emplace(funcs,                                                                             \
                                    "ER2007",                                                                          \
                                    [](Dune::MPIHelper::MPICommunicator global,                                        \
                                       Dune::MPIHelper::MPICommunicator local,                                         \
                                       Dune::XT::Common::Configuration config_in) {                                    \
                                      return new DMP::ER2007::FunctionName(global, local, config_in);                  \
                                    });                                                                                \
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
FunctionType* make_f(const std::map<std::string,
                                    std::function<FunctionType*(Dune::MPIHelper::MPICommunicator,
                                                                Dune::MPIHelper::MPICommunicator,
                                                                Dune::XT::Common::Configuration)>>& rets,
                     std::string name,
                     Dune::MPIHelper::MPICommunicator global,
                     Dune::MPIHelper::MPICommunicator local,
                     Dune::XT::Common::Configuration config_in)
{
  const auto it = rets.find(name);
  if (it == rets.end())
    DUNE_THROW(Dune::InvalidStateException, "no data for Problem. (toggle PROBLEM_NINE_ONLY?)");
  const auto& el = it->second;
  return el(global, local, config_in);
}

Problem::ProblemContainer::ProblemContainer(MPIHelper::MPICommunicator global,
                                            MPIHelper::MPICommunicator local,
                                            Dune::XT::Common::Configuration config_in)
  : config_(config_in)
  , name_(config_.get("problem.name", "Synthetic"))
  , data_(make_f(ModelProblemData_map, name_, global, local, config_in))
  , source_(make_f(Source_map, name_, global, local, config_in))
  , exact_solution_(make_f(ExactSolution_map, name_, global, local, config_in))
  , diffusion_(make_f(Diffusion_map, name_, global, local, config_in))
  , dirichlet_(make_f(DirichletData_map, name_, global, local, config_in))
  , neumann_(make_f(NeumannData_map, name_, global, local, config_in))
{
  getMutableModelData().problem_init(*this, global, local);
}

const CommonTraits::FunctionBaseType& DMP::ProblemContainer::getSource() const
{
  return *source_;
}

const CommonTraits::FunctionBaseType& DMP::ProblemContainer::getExactSolution() const
{
  return *exact_solution_;
}

const Problem::IModelProblemData& DMP::ProblemContainer::getModelData() const
{
  return *data_;
}

Problem::IModelProblemData& DMP::ProblemContainer::getMutableModelData()
{
  return *data_;
}

const Problem::DiffusionBase& DMP::ProblemContainer::getDiffusion() const
{
  return *diffusion_;
}

Problem::DiffusionBase& DMP::ProblemContainer::getMutableDiffusion()
{
  const auto& diffusion = getDiffusion();
  return const_cast<Problem::DiffusionBase&>(diffusion);
}

const Problem::DirichletDataBase& DMP::ProblemContainer::getDirichletData() const
{
  return *dirichlet_;
}

const Problem::NeumannDataBase& DMP::ProblemContainer::getNeumannData() const
{
  return *neumann_;
}

const std::string DMP::ProblemContainer::name() const
{
  return name_;
}

const Dune::XT::Common::Configuration& Problem::ProblemContainer::config() const
{
  return config_;
}

Dune::XT::Common::Configuration& Problem::ProblemContainer::config()
{
  return config_;
}
