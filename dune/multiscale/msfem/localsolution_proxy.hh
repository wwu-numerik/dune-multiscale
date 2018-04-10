#ifndef DUNE_MULTISCALE_MSFEM_LOCALSOLUTION_PROXY_HH
#define DUNE_MULTISCALE_MSFEM_LOCALSOLUTION_PROXY_HH

#include <dune/stuff/functions.hh>
#include <dune/multiscale/common/traits.hh>
#include <dune/multiscale/msfem/msfem_traits.hh>
#include <dune/multiscale/msfem/localproblems/localgridsearch.hh>

#include <dune/xt/common/parallel/threadstorage.hh>
#include <dune/xt/common/configuration.hh>
#include <unordered_map>

namespace Dune {
namespace Multiscale {

namespace Problem {
class ProblemContainer;
}
class ProxyGridview;

class LocalGridSearch;

/**
 * Fake DiscreteFunction that forwards localFunction calls to appropriate local_correction
 **/
class LocalsolutionProxy : public MsFEMTraits::LocalGridConstDiscreteFunctionType
{
  typedef LocalsolutionProxy ThisType;
  typedef MsFEMTraits::LocalGridConstDiscreteFunctionType BaseType;
  typedef CommonTraits::GridType::Traits::LeafIndexSet LeafIndexSetType;
  typedef typename BaseType::LocalfunctionType LocalFunctionType;

public:
  typedef std::unordered_map<typename LeafIndexSetType::IndexType,
                             std::unique_ptr<MsFEMTraits::LocalGridDiscreteFunctionType>>
      CorrectionsMapType;

  LocalsolutionProxy(CorrectionsMapType&& corrections,
                     const CommonTraits::SpaceType& coarseSpace,
                     const LocalGridList& gridlist,
                     const DMP::ProblemContainer& problem);

  std::unique_ptr<LocalFunctionType> local_function(const typename BaseType::EntityType& entity) const;

  void add(const CommonTraits::DiscreteFunctionType& coarse_func);

  LocalGridSearch& search();
  void visualize_parts(const XT::Common::Configuration& config) const;

  void visualize(const std::string&) const;

private:
  CorrectionsMapType corrections_;
  const CommonTraits::GridViewType view_;
  const LeafIndexSetType& index_set_;
  Dune::XT::Common::PerThreadValue<LocalGridSearch> search_;
};

} // namespace Multiscale {
} // namespace Dune {

#endif // DUNE_MULTISCALE_MSFEM_LOCALSOLUTION_PROXY_HH
