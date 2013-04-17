// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef MSFEM_TRAITS_HH
#define MSFEM_TRAITS_HH

#include <dune/multiscale/common/traits.hh>
#include <dune/multiscale/tools/errorestimation/MsFEM/msfem_elliptic_error_estimator.hh>
#include <dune/multiscale/msfem/msfem_solver.hh>
#include <dune/multiscale/tools/subgrid_io.hh>
#include <dune/multiscale/msfem/localproblems/subgrid-list.hh>
#include <dune/subgrid/subgrid.hh>
#include <dune/fem/io/file/dataoutput.hh>


namespace Dune {
namespace Multiscale {
namespace MsFEM {

template < class T >
class MacroMicroGridSpecifier;

//! type construction for the MSFEM code
struct MsFEMTraits {
  typedef MacroMicroGridSpecifier< typename CommonTraits::DiscreteFunctionSpaceType >                          MacroMicroGridSpecifierType;
  typedef Dune::SubGrid< CommonTraits::GridType::dimension, typename CommonTraits::GridType >                  SubGridType;
  typedef SubGridList< typename CommonTraits::DiscreteFunctionType, SubGridType, MacroMicroGridSpecifierType > SubGridListType;

  //! -------------------------- MsFEM error estimator ----------------------------
  typedef MsFEMErrorEstimator< typename CommonTraits::DiscreteFunctionType,
                               typename CommonTraits::DiffusionType,
                               typename CommonTraits::FirstSourceType,
                               MacroMicroGridSpecifierType,
                               SubGridListType >
    MsFEMErrorEstimatorType;
  //! -----------------------------------------------------------------------------
};

} //namespace MsFEM {
} //namespace Multiscale {
} //namespace Dune {}

#endif // MSFEM_TRAITS_HH
