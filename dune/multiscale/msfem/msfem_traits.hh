#ifndef MSFEM_TRAITS_HH
#define MSFEM_TRAITS_HH

#include <dune/multiscale/tools/errorestimation/MsFEM/msfem_elliptic_error_estimator.hh>
#include <dune/multiscale/tools/solver/MsFEM/msfem_solver.hh>
#include <dune/multiscale/tools/subgrid_io.hh>
#include <dune/multiscale/tools/solver/MsFEM/msfem_localproblems/subgrid-list.hh>
#include <dune/subgrid/subgrid.hh>
#include <dune/fem/io/file/dataoutput.hh>

namespace Dune {
namespace Multiscale {
namespace MsFEM {

template < class T >
class MacroMicroGridSpecifier;

//! --------- typedefs for the coefficient and data functions ------------------
struct MsfemTraits {
  //! ----- typedefs for the macro grid and the corresponding discrete space -----
  typedef Dune::GridSelector::GridType
  GridType;
  // Dune::InteriorBorder_Partition or Dune::All_Partition >?
  // see:
  // http://www.dune-project.org/doc/doxygen/dune-grid-html/group___g_i_related_types.html#ga5b9e8102d7f70f3f4178182629d98b6
  typedef Dune::AdaptiveLeafGridPart< GridType /*,Dune::All_Partition*/ > GridPartType;

  typedef Dune::GridPtr< GridType > GridPointerType;

  typedef Dune::FunctionSpace< double, double, WORLDDIM, 1 > FunctionSpaceType;
  // type of first source term (right hand side of differential equation or type of 'f')
  typedef Problem::FirstSource< FunctionSpaceType > FirstSourceType;

  // type of (possibly non-linear) diffusion term (i.e. 'A^{\epsilon}')
  typedef Problem::Diffusion< FunctionSpaceType > DiffusionType;

  // default type for any missing coefficient function (e.g. advection,...)
  typedef Problem::DefaultDummyFunction< FunctionSpaceType > DefaultDummyFunctionType;

  // type of exact solution (in general unknown)
  typedef Problem::ExactSolution< FunctionSpaceType > ExactSolutionType;
  typedef Dune::GridFunctionAdapter< ExactSolutionType, GridPartType >
    DiscreteExactSolutionType;     // for data output with paraview or grape


  //! ----  typedefs for the standard discrete function space (macroscopic) -----
  typedef FunctionSpaceType::DomainType DomainType;
  //! define the type of elements of the codomain v(\Omega) (typically a subset of \R)
  typedef FunctionSpaceType::RangeType RangeType;
  typedef std::vector< RangeType > RangeVector;
  typedef std::vector< RangeVector > RangeVectorVector;
  //! defines the function space to which the numerical solution belongs to
  //! see dune/fem/lagrangebase.hh
  typedef Dune::LagrangeDiscreteFunctionSpace< FunctionSpaceType, GridPartType, POLORDER >  // 1 = POLORDER
    DiscreteFunctionSpaceType;

  typedef Dune::AdaptiveDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;
  //!-----------------------------------------------------------------------------

  //!------------------------- for adaptive grid refinement ---------------------------------
  //! For adaption:
  //! type of restrict-prolong operator
  typedef Dune::RestrictProlongDefault< DiscreteFunctionType >
    RestrictProlongOperatorType;
  //! type of the adaption manager
  typedef Dune::AdaptationManager< GridType, RestrictProlongOperatorType >
    AdaptationManagerType;
  //!---------------------------------------------------------------------------------------

  typedef MacroMicroGridSpecifier< DiscreteFunctionSpaceType >                          MacroMicroGridSpecifierType;
  typedef Dune::SubGrid< GridType::dimension, GridType >                                      SubGridType;
  typedef Dune::SubGridList< DiscreteFunctionType, SubGridType, MacroMicroGridSpecifierType > SubGridListType;

  //! -------------------------- MsFEM error estimator ----------------------------
  typedef MsFEMErrorEstimator< DiscreteFunctionType,
                               DiffusionType,
                               FirstSourceType,
                               MacroMicroGridSpecifierType,
                               SubGridListType >
  MsFEMErrorEstimatorType;
  //! -----------------------------------------------------------------------------

  //! ------------------ typedefs and classes for data output ---------------------
  typedef Dune::tuple< const DiscreteFunctionType* >      IOTupleType;
  typedef Dune::DataOutput< GridType, IOTupleType > DataOutputType;
  // just for the discretized exact solution (in case it is available)
  typedef Dune::tuple< const DiscreteExactSolutionType* > ExSolIOTupleType;
  // just for the discretized exact solution (in case it is available)
  typedef Dune::DataOutput< GridType, ExSolIOTupleType > ExSolDataOutputType;
};

} //namespace MsFEM {
} //namespace Multiscale {
} //namespace Dune {}

#endif // MSFEM_TRAITS_HH
