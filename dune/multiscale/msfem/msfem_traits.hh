// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DHUUNE_MSFEM_TRAITS_HH
#define DHUUNE_MSFEM_TRAITS_HH

#include <dune/multiscale/common/la_backend.hh>
#include <dune/multiscale/common/traits.hh>

#include <dune/grid/spgrid.hh>
#include <dune/fem/space/lagrange.hh>

namespace Dune {

namespace Multiscale {
namespace MsFEM {
class MacroMicroGridSpecifier;
class LocalGridList;
template <class T, class R, class S, class G, class H>
class MsFEMErrorEstimator;

// ! type construction for the MSFEM code
struct MsFEMTraits {
  typedef MacroMicroGridSpecifier MacroMicroGridSpecifierType;
  typedef typename CommonTraits::DiscreteFunctionType::DiscreteFunctionSpaceType::FunctionSpaceType FunctionSpaceType;
  typedef Dune::SPGrid<double, CommonTraits::GridType::dimension> LocalGridType;
  typedef Fem::AdaptiveLeafGridPart<LocalGridType> LocalGridPartType;
  typedef Fem::LagrangeDiscreteFunctionSpace<FunctionSpaceType, LocalGridPartType, st_lagrangespace_order> LocalGridDiscreteFunctionSpaceType;

  typedef typename BackendChooser<LocalGridDiscreteFunctionSpaceType>::DiscreteFunctionType LocalGridDiscreteFunctionType;
  typedef LocalGridList LocalGridListType;

  // ! -------------------------- MsFEM error estimator ----------------------------
  typedef MsFEMErrorEstimator<typename CommonTraits::DiscreteFunctionType, typename CommonTraits::DiffusionType,
                              typename CommonTraits::FirstSourceType, MacroMicroGridSpecifierType,
                              LocalGridListType> MsFEMErrorEstimatorType;
  // ! -----------------------------------------------------------------------------

  // the following two may change if we intend to use different meshes on coarse and fine level
  typedef typename CommonTraits::GridType::Codim<0>::Entity CoarseEntityType;
  typedef typename CommonTraits::DiscreteFunctionSpaceType::BasisFunctionSetType CoarseBaseFunctionSetType;
};

} // namespace MsFEM {
} // namespace Multiscale {
} // namespace Dune {}

namespace DMM = Dune::Multiscale::MsFEM;

#endif // DUNE_MSFEM_TRAITS_HH
