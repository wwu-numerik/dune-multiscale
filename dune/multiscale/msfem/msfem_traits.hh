// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DHUUNE_MSFEM_TRAITS_HH
#define DHUUNE_MSFEM_TRAITS_HH

#include <dune/multiscale/common/la_backend.hh>
#include <dune/multiscale/common/traits.hh>

#include <dune/grid/spgrid.hh>
#include <dune/grid/sgrid.hh>
#include <dune/stuff/functions/constant.hh>
#include <vector>

namespace Dune {
namespace Multiscale {
namespace MsFEM {

//! type construction for the MSFEM code
struct MsFEMTraits {
  typedef Dune::SPGrid<double, CommonTraits::GridType::dimension, SPIsotropicRefinement, No_Comm> LocalGridType;
  // change dirichletconstraints.cc bottom accordingly
  // typedef Dune::SGrid<CommonTraits::GridType::dimension, CommonTraits::GridType::dimension> LocalGridType;

  typedef DSG::Providers::ConstDefault<LocalGridType> LocalGridProviderType;
  typedef GDT::Spaces::ContinuousLagrangeProvider<LocalGridType, DSG::ChooseLayer::leaf, CommonTraits::gdt_backend_type,
                                                  st_lagrangespace_order, CommonTraits::FieldType,
                                                  CommonTraits::dimRange> SpaceProviderType;

  //! \todo not correct
  typedef typename SpaceProviderType::Type LocalGridDiscreteFunctionSpaceType;
  typedef LocalGridDiscreteFunctionSpaceType LocalSpaceType;
  typedef typename LocalGridDiscreteFunctionSpaceType::EntityType LocalEntityType;

  typedef typename BackendChooser<LocalGridDiscreteFunctionSpaceType>::DiscreteFunctionType
  LocalGridDiscreteFunctionType;
  typedef typename BackendChooser<LocalGridDiscreteFunctionSpaceType>::ConstDiscreteFunctionType
  LocalGridConstDiscreteFunctionType;
  typedef Stuff::Functions::Constant<LocalEntityType, CommonTraits::FieldType, CommonTraits::dimDomain,
                                     CommonTraits::FieldType, CommonTraits::dimRange> LocalConstantFunctionType;
  typedef typename LocalGridDiscreteFunctionSpaceType::GridViewType LocalGridViewType;

  typedef typename CommonTraits::GridType::Codim<0>::Entity CoarseEntityType;
  typedef typename CommonTraits::DiscreteFunctionSpaceType::BaseFunctionSetType CoarseBaseFunctionSetType;

  typedef std::vector<std::shared_ptr<LocalGridDiscreteFunctionType>> LocalSolutionVectorType;
};

} // namespace MsFEM {
} // namespace Multiscale {
} // namespace Dune {}

namespace DMM = Dune::Multiscale::MsFEM;

#endif // DUNE_MSFEM_TRAITS_HH
