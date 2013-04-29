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

namespace FEM{
template <class R, class T, class S>
class DiscreteEllipticOperator;
}
namespace HMM {

template <class R, class T, class S>
class ErrorEstimator;

//! type construction for the HMM algorithm
struct HMMTraits {

  //! --------- typedefs for the periodic micro grid and the corresponding discrete space ----
  typedef Dune::PeriodicLeafGridPart< typename CommonTraits::GridType > PeriodicGridPartType;
  typedef Dune::LagrangeDiscreteFunctionSpace< typename CommonTraits::FunctionSpaceType, PeriodicGridPartType, 1 > // 1 =POLORDER
    PeriodicDiscreteFunctionSpaceType;
  typedef Dune::AdaptiveDiscreteFunction< PeriodicDiscreteFunctionSpaceType > PeriodicDiscreteFunctionType;

  /** \brief --------------- solver for the linear system of equations ----------------------------
     * use Bi CG Stab [OEMBICGSTABOp] or GMRES [OEMGMRESOp] for non-symmetric matrices and CG [CGInverseOp] for symmetric
     ****ones. GMRES seems to be more stable, but is extremely slow!
     */
  typedef Dune::
    OEMBICGSQOp
  //  OEMBICGSTABOp
  //    OEMGMRESOp
    < typename CommonTraits::DiscreteFunctionType, CommonTraits::FEMMatrix > InverseFEMMatrix;

  //! --------------- the discrete operators (standard FEM and HMM) ------------------------
  //! discrete elliptic operator (corresponds with FEM Matrix)
  typedef Dune::Multiscale::FEM::DiscreteEllipticOperator< typename CommonTraits::DiscreteFunctionType,
                                                           typename CommonTraits::DiffusionType,
                                                           typename CommonTraits::MassTermType >
    EllipticOperatorType;

  //! --------------------------------------------------------------------------------------


  //! --------------- ERROR ESTIMATOR NOT YET IMPLEMENTED ------------------------
    typedef ErrorEstimator< PeriodicDiscreteFunctionType,
                            typename CommonTraits::DiscreteFunctionType,
                            typename CommonTraits::DiffusionType > ErrorEstimatorType;
}; // struct  HMMTraits

} //namespace HMM {
} //namespace Multiscale {
} //namespace Dune {

#endif // DUNE_MS_HMM_TYPES_HH
