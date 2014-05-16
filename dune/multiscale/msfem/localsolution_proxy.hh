#ifndef DUNE_MULTISCALE_MSFEM_LOCALSOLUTION_PROXY_HH
#define DUNE_MULTISCALE_MSFEM_LOCALSOLUTION_PROXY_HH

#include <dune/stuff/functions.hh>
#include <dune/multiscale/common/traits.hh>

namespace Dune {
namespace Multiscale {
namespace MsFEM {

template <class T>
class LocalsolutionProxy;

/**
 * Fake DiscreteFunction that forwards localFunction calls to appropriate local_correction
 **/
template <class SearchType>
class LocalsolutionProxy : public MsFEMTraits::LocalGridConstDiscreteFunctionType
{
  typedef LocalsolutionProxy<SearchType> ThisType;
  typedef MsFEMTraits::LocalGridConstDiscreteFunctionType BaseType;
  typedef CommonTraits::GridType::Traits::LeafIndexSet LeafIndexSetType;
  typedef typename BaseType::LocalfunctionType LocalFunctionType;

public:
  typedef std::map<typename LeafIndexSetType::IndexType,
                   std::unique_ptr<DMM::MsFEMTraits::LocalGridDiscreteFunctionType>> CorrectionsMapType;

  LocalsolutionProxy(const CorrectionsMapType& corrections, const LeafIndexSetType& index_set,
                     SearchType& search)
    : BaseType(*corrections.begin()->second)
    , corrections_(corrections)
    , index_set_(index_set)
    , search_(search) {}

  std::unique_ptr<LocalFunctionType> local_function(const typename BaseType::EntityType& entity) const {
    const auto& coarse_cell = *search_.current_coarse_pointer();
    auto it = corrections_.find(index_set_.index(coarse_cell));
    if (it != corrections_.end())
      return it->second->local_function(entity);
    DUNE_THROW(InvalidStateException, "Coarse cell was not found!");
  }

  virtual ThisType* copy() const
  {
    DUNE_THROW(NotImplemented, "");
  }

private:
  const CorrectionsMapType& corrections_;
  const LeafIndexSetType& index_set_;
  SearchType& search_;
};

} // namespace MsFEM {
} // namespace Multiscale {
} // namespace Dune {

#endif // DUNE_MULTISCALE_MSFEM_LOCALSOLUTION_PROXY_HH
