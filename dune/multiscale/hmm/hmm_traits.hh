#ifndef DUNE_MS_HMM_TYPES_HH
#define DUNE_MS_HMM_TYPES_HH

#ifdef HAVE_CMAKE_CONFIG
 #include "cmake_config.h"
#elif defined (HAVE_CONFIG_H)
 #include <config.h>
#endif // ifdef HAVE_CMAKE_CONFIG

#include <dune/common/tuples.hh>

#include <dune/multiscale/tools/solver/HMM/cell_problem_solving/numbering.hh>
#include <dune/multiscale/tools/solver/HMM/cell_problem_solving/discreteoperator.hh>

#include <dune/multiscale/tools/assembler/matrix_assembler/elliptic_fem_matrix_assembler.hh>
#include <dune/multiscale/tools/assembler/matrix_assembler/elliptic_hmm_matrix_assembler.hh>
#include <dune/multiscale/tools/errorestimation/HMM/elliptic_error_estimator.hh>
#include <dune/multiscale/tools/homogenizer/elliptic_homogenizer.hh>

#include <dune/fem/operator/2order/lagrangematrixsetup.hh>
#include <dune/fem/space/common/adaptmanager.hh>
#include <dune/fem/io/file/dataoutput.hh>

struct HMMTraits {
//! --------- typedefs for the macro grid and the corresponding discrete space -------------
typedef Dune::GridSelector::GridType
    GridType;
// Dune::InteriorBorder_Partition or Dune::All_Partition >?
// see:
// http://www.dune-project.org/doc/doxygen/dune-grid-html/group___g_i_related_types.html#ga5b9e8102d7f70f3f4178182629d98b6
typedef Dune::AdaptiveLeafGridPart< GridType /*,Dune::All_Partition*/ > GridPartType;
typedef Dune::GridPtr< GridType >                                       GridPointerType;
typedef Dune::FunctionSpace< double, double, WORLDDIM, 1 >              FunctionSpaceType;
//!-----------------------------------------------------------------------------------------

//! --------- typedefs for the coefficient and data functions ------------------------------
typedef Problem::ModelProblemData ModelProblemDataType;
// type of first source term (right hand side of differential equation or type of 'f')
typedef Problem::FirstSource< FunctionSpaceType > FirstSourceType;
// type of second source term 'G' (second right hand side of differential equation 'div G')
typedef Problem::SecondSource< FunctionSpaceType > SecondSourceType;
// type of (possibly non-linear) diffusion term (i.e. 'A^{\epsilon}')
typedef Problem::Diffusion< FunctionSpaceType > DiffusionType;
// type of mass (or reaction) term (i.e. 'm' or 'c')
typedef Problem::MassTerm< FunctionSpaceType > MassTermType;
// default type for any missing coefficient function (e.g. advection,...)
typedef Problem::DefaultDummyFunction< FunctionSpaceType > DefaultDummyFunctionType;
typedef Dune::Homogenizer< GridType, DiffusionType > HomogenizerType;
typedef Problem::ConstantDiffusionMatrix< FunctionSpaceType, HomogenizerType::HomTensorType >
    HomDiffusionType;
//!-----------------------------------------------------------------------------------------

//! ---------  typedefs for the standard discrete function space (macroscopic) -------------

// type of exact solution (in general unknown)
typedef Problem::ExactSolution< FunctionSpaceType >            ExactSolutionType;
typedef Dune::GridFunctionAdapter< ExactSolutionType, GridPartType > DiscreteExactSolutionType;     // for data output with
                                                                                              // paraview or grape

typedef FunctionSpaceType::DomainType DomainType;
//! define the type of elements of the codomain v(\Omega) (typically a subset of \R)
typedef FunctionSpaceType::RangeType RangeType;
//! defines the function space to which the numerical solution belongs to
//! see dune/fem/lagrangebase.hh
typedef Dune::LagrangeDiscreteFunctionSpace< FunctionSpaceType, GridPartType, 1 >  // 1=POLORDER
    DiscreteFunctionSpaceType;
typedef DiscreteFunctionSpaceType::DomainFieldType            TimeType;
typedef DiscreteFunctionSpaceType::JacobianRangeType          JacobianRangeType;
typedef GridType::Codim< 0 >::Entity                          EntityType;
typedef GridType::Codim< 0 >::EntityPointer                   EntityPointerType;
typedef GridType::Codim< 0 >::Geometry                        EntityGeometryType;
typedef GridType::Codim< 1 >::Geometry                        FaceGeometryType;
typedef DiscreteFunctionSpaceType::BaseFunctionSetType        BaseFunctionSetType;
typedef Dune::CachingQuadrature< GridPartType, 0 >                  EntityQuadratureType;
typedef Dune::CachingQuadrature< GridPartType, 1 >                  FaceQuadratureType;
typedef DiscreteFunctionSpaceType::RangeFieldType             RangeFieldType;
typedef Dune::AdaptiveDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;
typedef DiscreteFunctionType::LocalFunctionType               LocalFunctionType;
typedef DiscreteFunctionType::DofIteratorType                 DofIteratorType;
//!-----------------------------------------------------------------------------------------

//! --------- typedefs for the periodic micro grid and the corresponding discrete space ----
typedef Dune::PeriodicLeafGridPart< GridType > PeriodicGridPartType;
typedef Dune::LagrangeDiscreteFunctionSpace< FunctionSpaceType, PeriodicGridPartType, 1 > // 1 =POLORDER
PeriodicDiscreteFunctionSpaceType;
typedef Dune::AdaptiveDiscreteFunction< PeriodicDiscreteFunctionSpaceType > PeriodicDiscreteFunctionType;
//!-----------------------------------------------------------------------------------------

//! ------------ cell problem solver and numbering manager -----------------------------------------
typedef Dune::CellProblemNumberingManager< DiscreteFunctionSpaceType > CellProblemNumberingManagerType;


//! --------------------- the standard matrix traits -------------------------------------
struct MatrixTraits
{
  typedef DiscreteFunctionSpaceType                          RowSpaceType;
  typedef DiscreteFunctionSpaceType                          ColumnSpaceType;
  typedef Dune::LagrangeMatrixSetup< false >                       StencilType;
  typedef Dune::ParallelScalarProduct< DiscreteFunctionSpaceType > ParallelScalarProductType;

  template< class M >
  struct Adapter
  {
    typedef Dune::LagrangeParallelMatrixAdapter< M > MatrixAdapterType;
  };
};

//! --------------------- type of fem stiffness matrix -----------------------------------
typedef Dune::SparseRowMatrixOperator< DiscreteFunctionType, DiscreteFunctionType, MatrixTraits > FEMMatrix;

/** \brief --------------- solver for the linear system of equations ----------------------------
   * use Bi CG Stab [OEMBICGSTABOp] or GMRES [OEMGMRESOp] for non-symmetric matrices and CG [CGInverseOp] for symmetric
   ****ones. GMRES seems to be more stable, but is extremely slow!
   */
typedef Dune::
  OEMBICGSQOp
//  OEMBICGSTABOp
//    OEMGMRESOp
  < DiscreteFunctionType, FEMMatrix > InverseFEMMatrix;

//! --------------- the discrete operators (standard FEM and HMM) ------------------------
//! discrete elliptic operator (corresponds with FEM Matrix)
typedef Dune::DiscreteEllipticOperator< DiscreteFunctionType, DiffusionType, MassTermType > EllipticOperatorType;
// discrete elliptic HMM operator (corresponds with HMM (or HMFEM) Matrix)
typedef Dune::DiscreteEllipticHMMOperator< DiscreteFunctionType, PeriodicDiscreteFunctionType, DiffusionType,
                                     CellProblemNumberingManagerType > EllipticHMMOperatorType;
//! --------------------------------------------------------------------------------------

//! --------------- ERROR ESTIMATOR NOT YET IMPLEMENTED ------------------------
typedef Dune::ErrorEstimator< PeriodicDiscreteFunctionType,
                        DiscreteFunctionType,
                        DiffusionType > ErrorEstimatorType;

//! --------- typedefs and classes for data output -----------------------------------------
typedef Dune::tuple< DiscreteFunctionType* >      IOTupleType;
typedef Dune::DataOutput< GridType, IOTupleType > DataOutputType;

// just for the discretized exact solution (in case it is available)
typedef Dune::tuple< DiscreteExactSolutionType* > ExSolIOTupleType;
// just for the discretized exact solution (in case it is available)
typedef Dune::DataOutput< GridType, ExSolIOTupleType > ExSolDataOutputType;

//!------------------------- for adaptive grid refinement ---------------------------------
//! type of restrict-prolong operator
typedef Dune::RestrictProlongDefault< DiscreteFunctionType > RestrictProlongOperatorType;
//! type of the adaption manager
typedef Dune::AdaptationManager< GridType, RestrictProlongOperatorType > AdaptationManagerType;
//!---------------------------------------------------------------------------------------

static const int assembler_order = 2* DiscreteFunctionSpaceType::polynomialOrder + 2;
}; // struct  HMMTraits


#endif // DUNE_MS_HMM_TYPES_HH
