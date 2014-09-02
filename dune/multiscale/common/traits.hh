// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_MULTISCALE_COMMON_TRAITS_HH
#define DUNE_MULTISCALE_COMMON_TRAITS_HH

#include <dune/multiscale/common/la_backend.hh>
#include <dune/common/tuples.hh>
#include <dune/grid/spgrid.hh>
#include <dune/gdt/spaces/continuouslagrange.hh>
#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/operators/elliptic-cg.hh>
#include <dune/stuff/la/container.hh>
#include <dune/stuff/aliases.hh>
#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/functions/constant.hh>
#include <dune/stuff/grid/provider.hh>

namespace Dune {

template <class T>
struct GridPtr;

namespace GDT {
template <class T, class R>
class DiscreteFunction;
template <class T, class R>
class ConstDiscreteFunction;
}
namespace Multiscale {
namespace Problem {

struct DiffusionBase;
struct LowerOrderBase;
class NeumannDataBase;
class DirichletDataBase;
class IModelProblemData;

} // namespace Problem

//! Common Types, duh
struct CommonTraits {

  static constexpr int dimRange = 1;
  static constexpr int dimDomain = GridSelector::dimgrid;
  static constexpr int world_dim = GridSelector::dimworld;
  static_assert(dimDomain == world_dim, "we really don't want to use an embedded grid");
  typedef Dune::SPGrid<double, world_dim> GridType;
  //    typedef Dune::SGrid<world_dim, world_dim> GridType;
  //    typedef Dune::YaspGrid<world_dim> GridType;

  static constexpr unsigned int exact_solution_space_order = 3 * st_lagrangespace_order;

  typedef GridType::Codim<0>::Entity EntityType;
  typedef double FieldType;

  typedef DSG::Providers::ConstDefault<GridType> GridProviderType;

  static constexpr auto gdt_backend_type =
#if DUNE_MULTISCALE_WITH_DUNE_FEM
      GDT::ChooseSpaceBackend::fem;
#else
      GDT::ChooseSpaceBackend::pdelab;
#endif

  typedef GDT::Spaces::ContinuousLagrangeProvider<GridType, DSG::ChooseLayer::leaf, gdt_backend_type,
                                                  st_lagrangespace_order, FieldType, dimRange> SpaceProviderType;

  static constexpr auto st_gdt_grid_level = 0;
  typedef SpaceProviderType::Type GdtSpaceType;
  typedef GdtSpaceType DiscreteFunctionSpaceType;
  typedef GdtSpaceType::GridViewType GridViewType;

  typedef BackendChooser<DiscreteFunctionSpaceType>::LinearOperatorType LinearOperatorType;
  typedef BackendChooser<DiscreteFunctionSpaceType>::GdtVectorType GdtVectorType;
  typedef BackendChooser<DiscreteFunctionSpaceType>::DiscreteFunctionDataType DiscreteFunctionDataType;
  typedef BackendChooser<DiscreteFunctionSpaceType>::DiscreteFunctionType DiscreteFunctionType;
  typedef BackendChooser<DiscreteFunctionSpaceType>::ConstDiscreteFunctionType ConstDiscreteFunctionType;

  typedef Stuff::GlobalFunctionInterface<EntityType, FieldType, dimDomain, FieldType, dimRange> FunctionBaseType;
  typedef Stuff::GlobalFunctionInterface<EntityType, FieldType, dimDomain, FieldType, dimDomain, dimDomain>
  DiffusionFunctionBaseType;
  typedef Stuff::Functions::Constant<EntityType, FieldType, dimDomain, FieldType, dimRange> ConstantFunctionBaseType;
  typedef ConstantFunctionBaseType GdtConstantFunctionType;

  typedef DiscreteFunctionSpaceType::BaseFunctionSetType BaseFunctionSetType;
  typedef DiscreteFunctionSpaceType::DomainType DomainType;
  typedef DiscreteFunctionSpaceType::BaseFunctionSetType::RangeType RangeType;
  typedef FieldType TimeType;
  typedef DiscreteFunctionSpaceType::BaseFunctionSetType::JacobianRangeType JacobianRangeType;

  typedef GridType::Codim<0>::EntityPointer EntityPointerType;
  typedef GridType::Codim<0>::Geometry EntityGeometryType;
  typedef GridType::Codim<1>::Geometry FaceGeometryType;

  typedef FieldType RangeFieldType;
  typedef FieldType DomainFieldType;

  typedef std::shared_ptr<DiscreteFunctionType> DiscreteFunction_ptr;

  typedef std::vector<RangeType> RangeVector;
  typedef std::vector<RangeVector> RangeVectorVector;

  static constexpr int polynomial_order = DiscreteFunctionSpaceType::polOrder;
  static constexpr int quadrature_order = 2 * polynomial_order + 2;
};

template <class T = CommonTraits::DiscreteFunctionType>
std::shared_ptr<T> make_df_ptr(const std::string name, const typename T::SpaceType& space) {
  return std::make_shared<T>(space, name);
  //  return DSC::make_unique<T>(name, space);
}

} // namespace Multiscale
} // namespace Dune

#endif // DUNE_MULTISCALE_COMMON_TRAITS_HH
