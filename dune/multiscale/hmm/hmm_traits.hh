// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_MS_HMM_TYPES_HH
#define DUNE_MS_HMM_TYPES_HH

#include <dune/multiscale/common/traits.hh>

#include <dune/fem/gridpart/periodicgridpart/periodicgridpart.hh>
#include <dune/fem/operator/2order/lagrangematrixsetup.hh>
#include <dune/fem/operator/matrix/spmatrix.hh>
#include <dune/fem/space/common/adaptmanager.hh>
#include <dune/fem/io/file/dataoutput.hh>

namespace Dune {
namespace Multiscale {

namespace FEM {
class DiscreteEllipticOperator;
}
namespace HMM {

template <class R, class T, class S>
class ErrorEstimator;

//! type construction for the HMM algorithm
struct HMMTraits {

  //! --------- typedefs for the periodic micro grid and the corresponding discrete space ----
  typedef Dune::Fem::PeriodicLeafGridPart<typename CommonTraits::GridType> PeriodicGridPartType;
  typedef Dune::Fem::LagrangeDiscreteFunctionSpace<typename CommonTraits::FunctionSpaceType, PeriodicGridPartType,
                                                   st_lagrangespace_order> PeriodicDiscreteFunctionSpaceType;

  typedef typename BackendChooser<PeriodicDiscreteFunctionSpaceType>::DiscreteFunctionType PeriodicDiscreteFunctionType;

  typedef typename BackendChooser<typename CommonTraits::DiscreteFunctionSpaceType>::InverseOperatorType
  InverseOperatorType;

  //! --------------- the discrete operators (standard FEM and HMM) ------------------------
  //! discrete elliptic operator (corresponds with FEM Matrix)
  typedef Dune::Multiscale::FEM::DiscreteEllipticOperator EllipticOperatorType;

  //! --------------------------------------------------------------------------------------

  //! --------------- ERROR ESTIMATOR NOT YET IMPLEMENTED ------------------------
  typedef ErrorEstimator<PeriodicDiscreteFunctionType, typename CommonTraits::DiscreteFunctionType,
                         typename CommonTraits::DiffusionType> ErrorEstimatorType;

  //! type of restrict-prolong operator
  typedef Dune::Fem::RestrictProlongDefault<typename CommonTraits::DiscreteFunctionType> RestrictProlongOperatorType;
  //! type of the adaption manager
  typedef Dune::Fem::AdaptationManager<typename CommonTraits::GridType, RestrictProlongOperatorType>
  AdaptationManagerType;
}; // struct  HMMTraits

} // namespace HMM {
} // namespace Multiscale {
} // namespace Dune {

#endif // DUNE_MS_HMM_TYPES_HH
