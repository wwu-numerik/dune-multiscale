#ifndef DUNE_MULTISCALE_MSFEM_COARSE_SCALE_ASSEMBLER_HH
#define DUNE_MULTISCALE_MSFEM_COARSE_SCALE_ASSEMBLER_HH

#include <ostream>
#include <type_traits>
#include <assert.h>
#include <boost/noncopyable.hpp>

#include <dune/common/fmatrix.hh>
#include <dune/stuff/common/logging.hh>
#include <dune/multiscale/common/traits.hh>
#include <dune/multiscale/msfem/msfem_traits.hh>

namespace Dune {
namespace Multiscale {

class MsFEMCodim0Integral;
class MsFemCodim0Matrix;
class LocalSolutionManager;
class LocalGridList;

class MsFEMCodim0IntegralTraits {
public:
  typedef MsFEMCodim0Integral derived_type;
};

class MsFemCodim0MatrixTraits {
public:
  typedef MsFEMCodim0Integral LocalOperatorType;
  typedef MsFemCodim0Matrix derived_type;
};

class MsFEMCodim0Integral // LocalOperatorType
    : public GDT::LocalOperator::Codim0Interface<MsFEMCodim0IntegralTraits> {
public:
  typedef MsFEMCodim0IntegralTraits Traits;

private:
  static const size_t numTmpObjectsRequired_ = 1;
  typedef Stuff::LocalfunctionSetInterface<CommonTraits::EntityType, CommonTraits::DomainFieldType,
                                           CommonTraits::dimDomain, CommonTraits::RangeFieldType,
                                           CommonTraits::dimRange, 1> AnsatzLocalfunctionSetInterfaceType;
  typedef AnsatzLocalfunctionSetInterfaceType TestLocalfunctionSetInterfaceType;

public:
  explicit MsFEMCodim0Integral(const size_t over_integrate = 0)
    : over_integrate_(over_integrate) {}

  size_t numTmpObjectsRequired() const;

  void apply(Multiscale::LocalSolutionManager& localSolutionManager,
             const MsFEMTraits::LocalEntityType& localGridEntity, const TestLocalfunctionSetInterfaceType& testBase,
             const AnsatzLocalfunctionSetInterfaceType& ansatzBase,
             Dune::DynamicMatrix<CommonTraits::RangeFieldType>& ret,
             std::vector<Dune::DynamicMatrix<CommonTraits::RangeFieldType>>& tmpLocalMatrices) const;

private:
  const size_t over_integrate_;
};

class MsFemCodim0Matrix // LocalAssemblerType
    {
public:
  typedef MsFemCodim0MatrixTraits Traits;
  typedef typename Traits::LocalOperatorType LocalOperatorType;

  MsFemCodim0Matrix(const LocalOperatorType& op, LocalGridList& localGridList);

  const LocalOperatorType& localOperator() const;

private:
  static const size_t numTmpObjectsRequired_ = 1;

public:
  std::vector<size_t> numTmpObjectsRequired() const;

  void
  assembleLocal(const CommonTraits::SpaceType& testSpace, const CommonTraits::SpaceType& ansatzSpace,
                const CommonTraits::EntityType& coarse_grid_entity, CommonTraits::LinearOperatorType& systemMatrix,
                std::vector<std::vector<Dune::DynamicMatrix<CommonTraits::RangeFieldType>>>& tmpLocalMatricesContainer,
                std::vector<Dune::DynamicVector<size_t>>& tmpIndicesContainer) const; // ... assembleLocal(...)

private:
  const LocalOperatorType& localOperator_;
  LocalGridList& localGridList_;
}; // class LocalAssemblerCodim0Matrix

} // namespace Multiscale {
} // namespace Dune {

#endif // DUNE_MULTISCALE_MSFEM_COARSE_SCALE_ASSEMBLER_HH
