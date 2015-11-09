// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DHUUNE_MSFEM_TRAITS_HH
#define DHUUNE_MSFEM_TRAITS_HH

#include <dune/multiscale/common/la_backend.hh>
#include <dune/multiscale/common/traits.hh>

#include <dune/grid/spgrid.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/stuff/functions/constant.hh>
#include <vector>

namespace Dune {
namespace Multiscale {

template <class GridType>
struct LocalGridChooser {
  typedef Dune::SPGrid<double, CommonTraits::GridType::dimension, SPIsotropicRefinement> type;
};

template <>
struct LocalGridChooser<
    Dune::YaspGrid<CommonTraits::world_dim, Dune::EquidistantOffsetCoordinates<double, CommonTraits::world_dim>>> {
  typedef Dune::YaspGrid<CommonTraits::world_dim, Dune::EquidistantOffsetCoordinates<double, CommonTraits::world_dim>>
      type;
};

//! type construction for the MSFEM code
struct MsFEMTraits {
  // change dirichletconstraints.cc bottom accordingly
  typedef typename LocalGridChooser<CommonTraits::GridType>::type LocalGridType;

  typedef SpaceChooser<LocalGridType, CommonTraits::FieldType, CommonTraits::dimRange> SpaceChooserType;
  typedef typename SpaceChooserType::Type LocalSpaceType;
  typedef typename LocalSpaceType::EntityType LocalEntityType;

  typedef typename BackendChooser<LocalSpaceType>::DiscreteFunctionType LocalGridDiscreteFunctionType;
  typedef typename BackendChooser<LocalSpaceType>::ConstDiscreteFunctionType LocalGridConstDiscreteFunctionType;
  typedef Stuff::Functions::Constant<LocalEntityType, CommonTraits::FieldType, CommonTraits::dimDomain,
                                     CommonTraits::FieldType, CommonTraits::dimRange> LocalConstantFunctionType;
  typedef typename LocalSpaceType::GridViewType LocalGridViewType;

  typedef typename CommonTraits::GridType::Codim<0>::Entity CoarseEntityType;
  typedef typename CommonTraits::SpaceType::BaseFunctionSetType CoarseBaseFunctionSetType;

  typedef std::vector<std::shared_ptr<LocalGridDiscreteFunctionType>> LocalSolutionVectorType;
};

} // namespace Multiscale {
} // namespace Dune {}

#endif // DUNE_MSFEM_TRAITS_HH
