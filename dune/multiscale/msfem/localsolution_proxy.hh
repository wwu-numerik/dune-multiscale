#ifndef DUNE_MULTISCALE_MSFEM_LOCALSOLUTION_PROXY_HH
#define DUNE_MULTISCALE_MSFEM_LOCALSOLUTION_PROXY_HH

//#include <dune/fem/function/common/dofblock.hh>
#include <dune/fem/function/adaptivefunction.hh>

namespace Dune {

namespace Fem {
template <class T, int c>
class DofBlockProxy;
template <class T>
class Envelope;
}
namespace Multiscale {
namespace MsFEM {

template <class T>
class LocalsolutionProxy;

//! Traits class for AdaptiveDiscreteFunction and AdaptiveLocalFunction
template <class SearchType>
struct LocalsolutionProxyTraits {
  typedef MsFEMTraits::LocalGridDiscreteFunctionSpaceType DiscreteFunctionSpaceType;
  typedef LocalsolutionProxy<SearchType> DiscreteFunctionType;
  typedef MsFEMTraits::LocalGridDiscreteFunctionType ProxiedTraits;
  typedef LocalsolutionProxyTraits<SearchType> Traits;
  typedef typename ProxiedTraits::LocalFunctionStorageType LocalFunctionStorageType;
  typedef typename ProxiedTraits::LocalFunctionType LocalFunctionType;
  typedef typename DiscreteFunctionSpaceType::DomainType DomainType;
  typedef typename DiscreteFunctionSpaceType::RangeType RangeType;
  typedef typename DiscreteFunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename DiscreteFunctionSpaceType::RangeFieldType RangeFieldType;
  typedef RangeFieldType DofType;
  typedef typename DiscreteFunctionSpaceType::JacobianRangeType JacobianRangeType;
  typedef typename DiscreteFunctionSpaceType::BlockMapperType BlockMapperType;
  typedef typename DiscreteFunctionSpaceType::GridType GridType;
  typedef typename ProxiedTraits::MapperType MapperType;
  typedef typename ProxiedTraits::DofStorageType DofStorageType;
  typedef typename ProxiedTraits::MutableDofStorageType MutableDofStorageType;
  typedef typename DofStorageType::DofIteratorType DofIteratorType;
  typedef typename DofStorageType::ConstDofIteratorType ConstDofIteratorType;
  static const int blockSize = DiscreteFunctionSpaceType::localBlockSize;
  typedef typename ProxiedTraits::DofBlockType DofBlockType;
  typedef typename ProxiedTraits::ConstDofBlockType ConstDofBlockType;
  typedef typename ProxiedTraits::DofBlockPtrType DofBlockPtrType;
  typedef typename ProxiedTraits::ConstDofBlockPtrType ConstDofBlockPtrType;
}; // end struct LocalsolutionProxyTraits

template <class SearchType>
class LocalsolutionProxy : public Dune::Fem::DiscreteFunctionInterface<LocalsolutionProxyTraits<SearchType>> {
  typedef LocalsolutionProxyTraits<SearchType> TraitsType;
  typedef Dune::Fem::DiscreteFunctionInterface<TraitsType> BaseType;
  typedef CommonTraits::GridType::Traits::LeafIndexSet LeafIndexSetType;

  typedef typename BaseType::LocalFunctionType LocalFunctionType;

public:
  typedef typename TraitsType::DofType DofType;
  typedef std::map<typename LeafIndexSetType::IndexType,
                   std::unique_ptr<DMM::MsFEMTraits::LocalGridDiscreteFunctionType>> CorrectionsMapType;

  LocalsolutionProxy(const CorrectionsMapType& corrections, const LeafIndexSetType& index_set,
                     const MsFEMTraits::LocalGridDiscreteFunctionSpaceType& /*space_in*/, SearchType& search)
    : corrections_(corrections)
    , index_set_(index_set)
    , search_(search) {}

  const LocalFunctionType localFunction(const typename BaseType::EntityType& entity) const {
    const auto& coarse_cell = *search_.current_coarse_pointer();
    auto it = corrections_.find(index_set_.index(coarse_cell));
    if (it != corrections_.end())
      return it->second->localFunction(entity);
    assert(false);
  }

  virtual bool read_xdr(const std::string) { return false; }
  virtual bool write_xdr(const std::string) const { return false; }
  virtual bool read_ascii(const std::string) { return false; }
  virtual bool write_ascii(const std::string) const { return false; }

private:
  const CorrectionsMapType& corrections_;
  const LeafIndexSetType& index_set_;
  SearchType& search_;
};

} // namespace MsFEM {
} // namespace Multiscale {
} // namespace Dune {

#endif // DUNE_MULTISCALE_MSFEM_LOCALSOLUTION_PROXY_HH
