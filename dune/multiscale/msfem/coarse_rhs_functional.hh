// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_MSFEM_RIGHT_HAND_SIDE_FUNCTIONAL_HH
#define DUNE_MSFEM_RIGHT_HAND_SIDE_FUNCTIONAL_HH

#include <dune/multiscale/common/traits.hh>
#include <dune/multiscale/msfem/msfem_traits.hh>
#include <dune/multiscale/problems/base.hh>

#include <dune/gdt/functionals/base.hh>
#include <dune/gdt/localfunctional/codim0.hh>
#include <dune/gdt/localfunctional/codim1.hh>
#include <dune/gdt/localevaluation/product.hh>
#include <dune/gdt/assembler/system.hh>

namespace Dune {
namespace Multiscale {

class LocalGridList;
class LocalSolutionManager;
class RhsCodim0Integral;
class CoarseRhsFunctional;
class RhsCodim0Vector;

class RhsCodim0IntegralTraits // LocalOperator
    {
public:
  typedef RhsCodim0Integral derived_type;
};

class RhsCodim0VectorTraits {
public:
  typedef RhsCodim0Vector derived_type;
}; // class LocalAssemblerCodim0MatrixTraits

class RhsCodim0Integral // LocalFunctionalType
    : public GDT::LocalOperator::Codim0Interface<RhsCodim0IntegralTraits> {
public:
  typedef RhsCodim0IntegralTraits Traits;

private:
  static constexpr size_t numTmpObjectsRequired_ = 1;
  typedef Stuff::LocalfunctionSetInterface<CommonTraits::EntityType, CommonTraits::DomainFieldType,
                                           CommonTraits::dimDomain, CommonTraits::RangeFieldType,
                                           CommonTraits::dimRange, 1> TestLocalfunctionSetInterfaceType;

public:
  explicit RhsCodim0Integral(const size_t over_integrate = 0)
    : over_integrate_(over_integrate) {}

  size_t numTmpObjectsRequired() const;

  void apply(MsFEMTraits::LocalGridDiscreteFunctionType& dirichletExtension,
             Multiscale::LocalSolutionManager& localSolutionManager,
             const MsFEMTraits::LocalEntityType& localGridEntity, const TestLocalfunctionSetInterfaceType& testBase,
             Dune::DynamicVector<CommonTraits::RangeFieldType>& ret,
             std::vector<Dune::DynamicVector<CommonTraits::RangeFieldType>>& tmpLocalVectors) const;

private:
  const size_t over_integrate_;
};

class RhsCodim0Vector // LocalAssemblerType
    {
public:
  typedef RhsCodim0VectorTraits Traits;

  RhsCodim0Vector(const RhsCodim0Integral& func, LocalGridList& localGridList)
    : localFunctional_(func)
    , localGridList_(localGridList) {}

  const RhsCodim0Integral& localFunctional() const { return localFunctional_; }

private:
  static constexpr size_t numTmpObjectsRequired_ = 1;

public:
  std::vector<size_t> numTmpObjectsRequired() const;

  void
  assembleLocal(const CommonTraits::SpaceType& testSpace, const CommonTraits::EntityType& coarse_grid_entity,
                CommonTraits::GdtVectorType& systemVector,
                std::vector<std::vector<Dune::DynamicVector<CommonTraits::RangeFieldType>>>& tmpLocalVectorContainer,
                Dune::DynamicVector<size_t>& tmpIndices) const; // ... assembleLocal(...)

private:
  const RhsCodim0Integral& localFunctional_;
  LocalGridList& localGridList_;
}; // class RhsCodim0Vector

class CoarseRhsFunctionalTraits {
  typedef CommonTraits::GdtVectorType VectorImp;
  typedef CommonTraits::SpaceType SpaceImp;
  typedef typename CommonTraits::InteriorGridViewType GridViewImp;

public:
  typedef Problem::DiffusionBase FunctionType;

  static_assert(
      std::is_base_of<
          Stuff::LocalizableFunctionInterface<typename FunctionType::EntityType, typename FunctionType::DomainFieldType,
                                              FunctionType::dimDomain, typename FunctionType::RangeFieldType,
                                              FunctionType::dimRange, FunctionType::dimRangeCols>,
          FunctionType>::value,
      "FunctionType has to be derived from Stuff::LocalizableFunctionInterface!");
  static_assert(GDT::is_space< SpaceImp >::value,
                "SpaceImp has to be derived from SpaceInterface!");

  typedef GDT::SystemAssembler<SpaceImp, GridViewImp, SpaceImp> SystemAssemblerType;
  typedef CoarseRhsFunctional derived_type;
  typedef VectorImp VectorType;
  typedef SpaceImp SpaceType;
  typedef GridViewImp GridViewType;
  typedef typename VectorType::ScalarType ScalarType;
}; // class CoarseRhsFunctionalTraits

class CoarseRhsFunctional : public GDT::Functionals::VectorBased<CoarseRhsFunctionalTraits>,
                            public CoarseRhsFunctionalTraits::SystemAssemblerType {
  typedef GDT::Functionals::VectorBased<CoarseRhsFunctionalTraits> FunctionalBaseType;
  typedef CoarseRhsFunctionalTraits::SystemAssemblerType AssemblerBaseType;

  typedef RhsCodim0Integral LocalFunctionalType;
  typedef RhsCodim0Vector LocalAssemblerType;

public:
  typedef CoarseRhsFunctionalTraits Traits;
  typedef typename Traits::VectorType VectorType;
  typedef typename Traits::SpaceType SpaceType;
  typedef typename Traits::GridViewType GridViewType;

  CoarseRhsFunctional(VectorType& vec, const SpaceType& spc, LocalGridList& localGridList, const CommonTraits::InteriorGridViewType &interior);

  virtual ~CoarseRhsFunctional(){}

  virtual void assemble() override final;

private:
  const LocalFunctionalType local_functional_;
  const LocalAssemblerType local_assembler_;
}; // class CoarseRhsFunctional

} // end namespace Multiscale
} // end namespace Dune

#endif // DUNE_MSFEM_RIGHT_HAND_SIDE_FUNCTIONAL_HH
