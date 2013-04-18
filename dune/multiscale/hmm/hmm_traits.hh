// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_MS_HMM_TYPES_HH
#define DUNE_MS_HMM_TYPES_HH

#include <dune/multiscale/common/traits.hh>

#include <dune/fem/gridpart/periodicgridpart/periodicgridpart.hh>
#include <dune/multiscale/hmm/cell_problem_numbering.hh>
#include <dune/multiscale/hmm/discrete_cell_operator.hh>
#include <dune/multiscale/fem/elliptic_fem_matrix_assembler.hh>
#include <dune/multiscale/hmm/elliptic_hmm_matrix_assembler.hh>
#include <dune/multiscale/tools/errorestimation/HMM/elliptic_error_estimator.hh>

#include <dune/fem/operator/2order/lagrangematrixsetup.hh>
#include <dune/fem/space/common/adaptmanager.hh>
#include <dune/fem/io/file/dataoutput.hh>

namespace Dune {
namespace Multiscale {
namespace HMM {

//! type construction for the HMM algorithm
struct HMMTraits {

  //! --------- typedefs for the periodic micro grid and the corresponding discrete space ----
  typedef Dune::PeriodicLeafGridPart< typename CommonTraits::GridType > PeriodicGridPartType;
  typedef Dune::LagrangeDiscreteFunctionSpace< typename CommonTraits::FunctionSpaceType, PeriodicGridPartType, 1 > // 1 =POLORDER
    PeriodicDiscreteFunctionSpaceType;
  typedef Dune::AdaptiveDiscreteFunction< PeriodicDiscreteFunctionSpaceType > PeriodicDiscreteFunctionType;
  //!-----------------------------------------------------------------------------------------

  //! ------------ cell problem solver and numbering manager -----------------------------------------
  typedef Dune::CellProblemNumberingManager< typename CommonTraits::DiscreteFunctionSpaceType > CellProblemNumberingManagerType;


  //! --------------------- type of fem stiffness matrix -----------------------------------
  typedef Dune::SparseRowMatrixOperator< typename CommonTraits::DiscreteFunctionType, typename CommonTraits::DiscreteFunctionType, typename CommonTraits::MatrixTraits > FEMMatrix;

  /** \brief --------------- solver for the linear system of equations ----------------------------
     * use Bi CG Stab [OEMBICGSTABOp] or GMRES [OEMGMRESOp] for non-symmetric matrices and CG [CGInverseOp] for symmetric
     ****ones. GMRES seems to be more stable, but is extremely slow!
     */
  typedef Dune::
    OEMBICGSQOp
  //  OEMBICGSTABOp
  //    OEMGMRESOp
    < typename CommonTraits::DiscreteFunctionType, FEMMatrix > InverseFEMMatrix;

  //! --------------- the discrete operators (standard FEM and HMM) ------------------------
  //! discrete elliptic operator (corresponds with FEM Matrix)
  typedef Dune::Multiscale::FEM::DiscreteEllipticOperator< typename CommonTraits::DiscreteFunctionType,
                                                           typename CommonTraits::DiffusionType,
                                                           typename CommonTraits::MassTermType >
    EllipticOperatorType;
  // discrete elliptic HMM operator (corresponds with HMM (or HMFEM) Matrix)
  typedef Dune::DiscreteEllipticHMMOperator< typename CommonTraits::DiscreteFunctionType, PeriodicDiscreteFunctionType, typename CommonTraits::DiffusionType,
                                       CellProblemNumberingManagerType > EllipticHMMOperatorType;
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
