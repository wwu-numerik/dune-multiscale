#ifndef DUNE_MULTISCALE_MSFEM_LOCALSOLUTION_PROXY_HH
#define DUNE_MULTISCALE_MSFEM_LOCALSOLUTION_PROXY_HH

#include <dune/stuff/functions.hh>
#include <dune/multiscale/common/traits.hh>
#include <dune/multiscale/msfem/msfem_traits.hh>

namespace Dune {
namespace Multiscale {
namespace MsFEM {

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
  typedef std::map<typename LeafIndexSetType::IndexType,
                   std::unique_ptr<DMM::MsFEMTraits::LocalGridDiscreteFunctionType>> CorrectionsMapType;

  LocalsolutionProxy(const CorrectionsMapType& corrections, const LeafIndexSetType& index_set,
                     LocalGridSearch& search);

  std::unique_ptr<LocalFunctionType> local_function(const typename BaseType::EntityType& entity) const;

  virtual ThisType* copy() const
  {
    DUNE_THROW(NotImplemented, "");
  }

private:
  const CorrectionsMapType& corrections_;
  const LeafIndexSetType& index_set_;
  LocalGridSearch& search_;
};

} // namespace MsFEM {
} // namespace Multiscale {
} // namespace Dune {

#endif // DUNE_MULTISCALE_MSFEM_LOCALSOLUTION_PROXY_HH
