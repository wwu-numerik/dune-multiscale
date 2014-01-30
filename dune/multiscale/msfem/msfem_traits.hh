// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DHUUNE_MSFEM_TRAITS_HH
#define DHUUNE_MSFEM_TRAITS_HH

#include <dune/multiscale/common/la_backend.hh>
#include <dune/multiscale/common/traits.hh>

#include <dune/grid/spgrid.hh>
#include <dune/grid/sgrid.hh>
#include <dune/fem/space/lagrange.hh>
#include <vector>

namespace Dune {

namespace Multiscale {
namespace MsFEM {

//! type construction for the MSFEM code
struct MsFEMTraits {
  typedef typename CommonTraits::DiscreteFunctionType::DiscreteFunctionSpaceType::FunctionSpaceType FunctionSpaceType;
//  typedef Dune::SPGrid<double, CommonTraits::GridType::dimension> LocalGridType;
  typedef Dune::SGrid<CommonTraits::GridType::dimension, CommonTraits::GridType::dimension> LocalGridType;
  typedef Fem::AdaptiveLeafGridPart<LocalGridType> LocalGridPartType;
  typedef Fem::LagrangeDiscreteFunctionSpace<FunctionSpaceType, LocalGridPartType, st_lagrangespace_order> LocalGridDiscreteFunctionSpaceType;

  typedef typename LocalGridDiscreteFunctionSpaceType::IteratorType::Entity LocalEntityType;

  typedef typename BackendChooser<LocalGridDiscreteFunctionSpaceType>::DiscreteFunctionType LocalGridDiscreteFunctionType;

  typedef typename CommonTraits::GridType::Codim<0>::Entity CoarseEntityType;
  typedef typename CommonTraits::DiscreteFunctionSpaceType::BasisFunctionSetType CoarseBaseFunctionSetType;

  typedef std::vector<std::shared_ptr<LocalGridDiscreteFunctionType>> LocalSolutionVectorType;
};

} // namespace MsFEM {
} // namespace Multiscale {
} // namespace Dune {}

namespace DMM = Dune::Multiscale::MsFEM;

#endif // DUNE_MSFEM_TRAITS_HH
