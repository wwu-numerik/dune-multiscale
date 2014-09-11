#ifndef DUNE_MULTISCALE_MSFEM_LOCALSOLUTION_PROXY_HH
#define DUNE_MULTISCALE_MSFEM_LOCALSOLUTION_PROXY_HH

#include <dune/stuff/functions.hh>
#include <dune/multiscale/common/traits.hh>
#include <dune/multiscale/msfem/msfem_traits.hh>
#include <dune/multiscale/msfem/localproblems/localgridsearch.hh>
#include <unordered_map>

namespace Dune {
namespace Multiscale {

class ProxyGridview;

namespace MsFEM {

//class LocalGridSearch;
class LocalGridList;

/**
 * Fake DiscreteFunction that forwards localFunction calls to appropriate local_correction
 **/
class LocalsolutionProxy : public MsFEMTraits::LocalGridConstDiscreteFunctionType {
  typedef LocalsolutionProxy ThisType;
  typedef MsFEMTraits::LocalGridConstDiscreteFunctionType BaseType;
  typedef CommonTraits::GridType::Traits::LeafIndexSet LeafIndexSetType;
  typedef typename BaseType::LocalfunctionType LocalFunctionType;

public:
  typedef std::unordered_map<typename LeafIndexSetType::IndexType,
                             std::unique_ptr<DMM::MsFEMTraits::LocalGridDiscreteFunctionType>> CorrectionsMapType;

  LocalsolutionProxy(CorrectionsMapType &&corrections, const CommonTraits::DiscreteFunctionSpaceType &coarseSpace,
                     const LocalGridList& gridlist);

  std::unique_ptr<LocalFunctionType> local_function(const typename BaseType::EntityType& entity) const;

  ProxyGridview grid_view() const;

  LocalGridSearch& search() const;

private:
  const CorrectionsMapType corrections_;
  const LeafIndexSetType& index_set_;
  const std::unique_ptr<LocalGridSearch> search_;
  const LocalGridList& gridlist_;
};

} // namespace MsFEM {
} // namespace Multiscale {
} // namespace Dune {

#endif // DUNE_MULTISCALE_MSFEM_LOCALSOLUTION_PROXY_HH
