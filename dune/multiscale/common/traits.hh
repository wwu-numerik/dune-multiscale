// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_MULTISCALE_COMMON_TRAITS_HH
#define DUNE_MULTISCALE_COMMON_TRAITS_HH

#include <dune/common/tuples.hh>
#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/operators/elliptic-cg.hh>
#include <dune/gdt/spaces/cg.hh>
#include <dune/grid/spgrid.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/multiscale/common/la_backend.hh>
#include <dune/stuff/aliases.hh>
#include <dune/stuff/functions/constant.hh>
#include <dune/stuff/functions/expression.hh>
#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/grid/provider.hh>
#include <dune/stuff/la/container.hh>

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

template <class G, class R, int r>
struct SpaceChooser {
  static constexpr auto backend_type =
#if DUNE_MULTISCALE_WITH_DUNE_FEM
      GDT::ChooseSpaceBackend::fem;
#else
      GDT::ChooseSpaceBackend::pdelab;
#endif
  static constexpr auto partview_chooser = GDT::ChooseGridPartView<backend_type>::type;
  typedef DSG::LeafPartView<G, partview_chooser> PartViewType;
  typedef typename Stuff::Grid::Layer<G, DSG::ChooseLayer::leaf, partview_chooser>::Type GridLayerType;

private:
  typedef GDT::Spaces::CG::PdelabBased<GridLayerType, st_lagrangespace_order, R, r> PdelabType;
  typedef GDT::Spaces::CG::FemBased<GridLayerType, st_lagrangespace_order, R, r> FemType;

public:
  typedef typename std::conditional<(backend_type == GDT::ChooseSpaceBackend::fem), FemType, PdelabType>::type Type;
  static Type make_space(G& g) { return Type(PartViewType::create(g, 0)); }
  static Type make_space(GridLayerType& p) { return Type(p); }
};

//! Common Types, duh
struct CommonTraits {

  static constexpr int dimRange = 1;
  static constexpr int dimDomain = Dune::GridSelector::dimgrid;
  static constexpr int world_dim = Dune::GridSelector::dimworld;
  static_assert(dimDomain == world_dim, "we really don't want to use an embedded grid");
  typedef Dune::GridSelector::GridType GridType;

  static constexpr unsigned int exact_solution_space_order = 3 * st_lagrangespace_order;

  typedef GridType::Codim<0>::Entity EntityType;
  typedef double FieldType;

  typedef DSG::Providers::Default<GridType> GridProviderType;

  typedef SpaceChooser<GridType, FieldType, dimRange> SpaceChooserType;
  typedef typename SpaceChooserType::Type SpaceType;

  static constexpr auto st_gdt_grid_level = 0;

  typedef SpaceType::GridViewType GridViewType;
  typedef typename GridType::Partition<InteriorBorder_Partition>::LeafGridView InteriorGridViewType;
  static constexpr auto InteriorPartition = PartitionIteratorType::InteriorBorder_Partition;

  typedef BackendChooser<SpaceType>::LinearOperatorType LinearOperatorType;
  typedef BackendChooser<SpaceType>::GdtVectorType GdtVectorType;
  typedef BackendChooser<SpaceType>::DiscreteFunctionDataType DiscreteFunctionDataType;
  typedef BackendChooser<SpaceType>::DiscreteFunctionType DiscreteFunctionType;
  typedef BackendChooser<SpaceType>::ConstDiscreteFunctionType ConstDiscreteFunctionType;

  typedef Stuff::GlobalFunctionInterface<EntityType, FieldType, dimDomain, FieldType, dimRange> FunctionBaseType;
  typedef Stuff::GlobalFunctionInterface<EntityType, FieldType, dimDomain, FieldType, dimDomain, dimDomain>
      DiffusionFunctionBaseType;
  typedef Stuff::Functions::Constant<EntityType, FieldType, dimDomain, FieldType, dimRange> ConstantFunctionBaseType;
  typedef Stuff::Functions::Expression<EntityType, FieldType, dimDomain, FieldType, dimRange>
      ExpressionFunctionBaseType;
  typedef ConstantFunctionBaseType GdtConstantFunctionType;

  typedef SpaceType::BaseFunctionSetType BaseFunctionSetType;
  typedef SpaceType::DomainType DomainType;
  typedef SpaceType::BaseFunctionSetType::RangeType RangeType;
  typedef FieldType TimeType;
  typedef SpaceType::BaseFunctionSetType::JacobianRangeType JacobianRangeType;

  typedef GridType::Codim<0>::EntityPointer EntityPointerType;
  typedef GridType::Codim<0>::Geometry EntityGeometryType;
  typedef GridType::Codim<1>::Geometry FaceGeometryType;

  typedef FieldType RangeFieldType;
  typedef FieldType DomainFieldType;

  typedef std::shared_ptr<DiscreteFunctionType> DiscreteFunction_ptr;

  typedef std::vector<RangeType> RangeVector;
  typedef std::vector<RangeVector> RangeVectorVector;

  static constexpr int polynomial_order = SpaceType::polOrder;
  static constexpr int quadrature_order = 2 * polynomial_order + 2;
};

template <class T = CommonTraits::DiscreteFunctionType>
std::shared_ptr<T> make_df_ptr(const std::string name, const typename T::SpaceType& space) {
  return std::make_shared<T>(space, name);
  //  return Dune::XT::Common::make_unique<T>(name, space);
}

} // namespace Multiscale
} // namespace Dune

namespace DMP = Dune::Multiscale::Problem;

#endif // DUNE_MULTISCALE_COMMON_TRAITS_HH
