#ifndef DUNE_FEM_TYPES_HH
#define DUNE_FEM_TYPES_HH

#ifdef HAVE_CMAKE_CONFIG
 #include "cmake_config.h"
#elif defined (HAVE_CONFIG_H)
 #include <config.h>
#endif // ifdef HAVE_CMAKE_CONFIG

#include <dune/common/tuples.hh>

#include <dune/multiscale/tools/assembler/matrix_assembler/elliptic_fem_matrix_assembler.hh>
#include <dune/multiscale/tools/assembler/righthandside_assembler.hh>
#include <dune/fem/operator/2order/lagrangematrixsetup.hh>
#include <dune/fem/space/common/adaptmanager.hh>
#include <dune/fem/io/file/dataoutput.hh>

struct FEMTraits {
//! --------- typedefs for the grid and the corresponding discrete space -------------
typedef Dune::GridSelector::GridType
    GridType;
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
//!-----------------------------------------------------------------------------------------

//! ---------------  typedefs for the standard discrete function space ---------------------

// type of exact solution (in general unknown)
typedef Problem::ExactSolution< FunctionSpaceType >            ExactSolutionType;
typedef Dune::GridFunctionAdapter< ExactSolutionType, GridPartType > DiscreteExactSolutionType;     // for data output with


typedef FunctionSpaceType::DomainType DomainType;
//! define the type of elements of the codomain v(\Omega) (typically a subset of \R)
typedef FunctionSpaceType::RangeType RangeType;
//! defines the function space to which the numerical solution belongs to
//! see dune/fem/lagrangebase.hh
typedef Dune::LagrangeDiscreteFunctionSpace< FunctionSpaceType, GridPartType, 1 >  // 1=POLORDER
    DiscreteFunctionSpaceType;

//typedef DiscreteFunctionSpaceType::JacobianRangeType          JacobianRangeType;
//typedef DiscreteFunctionSpaceType::RangeFieldType             RangeFieldType;
typedef Dune::AdaptiveDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;
//!-----------------------------------------------------------------------------------------

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
typedef //Dune::OEMBICGSQOp
  Dune::OEMBICGSTABOp
//  OEMGMRESOp
  < DiscreteFunctionType, FEMMatrix > InverseFEMMatrix;

//! --------------- the discrete operators (standard FEM) ----------------------------------
//! discrete elliptic operator (corresponds with FEM Matrix)
typedef Dune::DiscreteEllipticOperator< DiscreteFunctionType, DiffusionType, MassTermType >
     EllipticOperatorType;
//! ----------------------------------------------------------------------------------------

//! --------- typedefs and classes for data output -----------------------------------------
typedef Dune::tuple< DiscreteFunctionType* >      IOTupleType;
typedef Dune::DataOutput< GridType, IOTupleType > DataOutputType;

// just for the discretized exact solution (in case it is available)
typedef Dune::tuple< DiscreteExactSolutionType* > ExSolIOTupleType;
// just for the discretized exact solution (in case it is available)
typedef Dune::DataOutput< GridType, ExSolIOTupleType > ExSolDataOutputType;

static const int assembler_order = 2* DiscreteFunctionSpaceType::polynomialOrder + 2;
}; // struct  FEMTraits


#endif // DUNE_FEM_TYPES_HH
