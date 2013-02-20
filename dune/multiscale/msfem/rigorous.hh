#ifndef DUNE_MULTISCALE_RIGOROUS_HH
#define DUNE_MULTISCALE_RIGOROUS_HH

#include <dune/fem/gridpart/common/gridpart.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/space/lagrangespace.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/operator/lagrangeinterpolation.hh>

// to display data with ParaView:
#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/fem/io/file/dataoutput.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/io/file/datawriter.hh>
#include <dune/fem/space/common/adaptmanager.hh>
#include <dune/fem/misc/l2error.hh>
#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/misc/h1norm.hh>

#include <dune/multiscale/problems/elliptic_problems/selector.hh>
#include <dune/multiscale/tools/solver/FEM/fem_solver.hh>
#include <dune/multiscale/tools/solver/MsFEM/msfem_localproblems/subgrid-list.hh>
#include <dune/multiscale/tools/solver/MsFEM/rigorous_msfem_solver.hh>
#include <dune/multiscale/tools/meanvalue.hh>
#include <dune/multiscale/tools/improved_l2error.hh>
#include <dune/multiscale/tools/misc/h1error.hh>
#include <dune/multiscale/tools/disc_func_writer/discretefunctionwriter.hh>
#include <dune/multiscale/tools/misc/outputparameter.hh>
#include <dune/multiscale/problems/elliptic_problems/selector.hh>
#include <dune/multiscale/msfem/msfem_traits.hh>

namespace Dune {
namespace Multiscale {
namespace MsFEM {

void solution_output(const MsFEMTraits::DiscreteFunctionType& msfem_solution,
                     const MsFEMTraits::DiscreteFunctionType& coarse_part_msfem_solution,
                     const MsFEMTraits::DiscreteFunctionType& fine_part_msfem_solution,
                     Dune::myDataOutputParameters& outputparam,
                     int& total_refinement_level_,
                     int& coarse_grid_level_)
{
  using namespace Dune;
  //! ----------------- writing data output MsFEM Solution -----------------
  // --------- VTK data output for MsFEM solution --------------------------
  // create and initialize output class
  MsFEMTraits::IOTupleType msfem_solution_series(&msfem_solution);
  const auto& gridPart = msfem_solution.space().gridPart();
  std::string outstring;
  outputparam.set_prefix("/msfem_solution");
  outstring = "msfem_solution";

  MsFEMTraits::DataOutputType msfem_dataoutput(gridPart.grid(), msfem_solution_series, outputparam);
  msfem_dataoutput.writeData( 1.0 /*dummy*/, outstring );

  // create and initialize output class
  MsFEMTraits::IOTupleType coarse_msfem_solution_series(&coarse_part_msfem_solution);

  outputparam.set_prefix("/coarse_part_msfem_solution");
  outstring = "coarse_part_msfem_solution";

  MsFEMTraits::DataOutputType coarse_msfem_dataoutput(gridPart.grid(), coarse_msfem_solution_series, outputparam);
  coarse_msfem_dataoutput.writeData( 1.0 /*dummy*/, outstring );

  // create and initialize output class
  MsFEMTraits::IOTupleType fine_msfem_solution_series(&fine_part_msfem_solution);

  outputparam.set_prefix("/fine_part_msfem_solution");
  // write data
  outstring = "fine_msfem_solution";

  MsFEMTraits::DataOutputType fine_msfem_dataoutput(gridPart.grid(), fine_msfem_solution_series, outputparam);
  fine_msfem_dataoutput.writeData( 1.0 /*dummy*/, outstring);

  // ----------------------------------------------------------------------
  // ---------------------- write discrete msfem solution to file ---------
  const std::string location = (boost::format("msfem_solution_discFunc_refLevel_%d_%d")
                                %  total_refinement_level_ % coarse_grid_level_).str();
  DiscreteFunctionWriter(location).append(msfem_solution);
  //! --------------------------------------------------------------------
}

void data_output(const MsFEMTraits::GridPartType& gridPart,
                 const MsFEMTraits::DiscreteFunctionSpaceType& discreteFunctionSpace_coarse,
                 Dune::myDataOutputParameters& outputparam )
{
  using namespace Dune;
  // sequence stamp
  //! --------------------------------------------------------------------------------------

  //! -------------------------- writing data output Exact Solution ------------------------
  if (Problem::ModelProblemData::has_exact_solution)
  {
    const MsFEMTraits::ExactSolutionType u;
    const MsFEMTraits::DiscreteExactSolutionType discrete_exact_solution("discrete exact solution ", u, gridPart);
    // create and initialize output class
    MsFEMTraits::ExSolIOTupleType exact_solution_series(&discrete_exact_solution);
    outputparam.set_prefix("/exact_solution");
    MsFEMTraits::ExSolDataOutputType exactsol_dataoutput(gridPart.grid(), exact_solution_series, outputparam);
    // write data
    exactsol_dataoutput.writeData( 1.0 /*dummy*/, "exact-solution" );
    // -------------------------------------------------------
  }
  //! --------------------------------------------------------------------------------------

  //! --------------- writing data output for the coarse grid visualization ------------------
  MsFEMTraits::DiscreteFunctionType coarse_grid_visualization("Visualization of the coarse grid",
                                                 discreteFunctionSpace_coarse);
  coarse_grid_visualization.clear();
  // -------------------------- data output -------------------------
  // create and initialize output class
  MsFEMTraits::IOTupleType coarse_grid_series(&coarse_grid_visualization);

  const auto coarse_grid_fname = (boost::format("/coarse_grid_visualization_")).str();
  outputparam.set_prefix(coarse_grid_fname);
  MsFEMTraits::DataOutputType coarse_grid_dataoutput(discreteFunctionSpace_coarse.gridPart().grid(), coarse_grid_series, outputparam);
  // write data
  coarse_grid_dataoutput.writeData( 1.0 /*dummy*/, coarse_grid_fname );
  // -------------------------------------------------------
  //! --------------------------------------------------------------------------------------
}


//! algorithm
void algorithm(const std::string& macroGridName,
               int& total_refinement_level_,
               int& coarse_grid_level_,
               int& number_of_layers_ ) {
  using namespace Dune;

  DSC_LOG_INFO << "loading dgf: " << macroGridName << std::endl;
  // we might use further grid parameters (depending on the grid type, e.g. Alberta), here we switch to default values
  // for the parameters:
  // create a grid pointer for the DGF file belongig to the macro grid:
  MsFEMTraits::GridPointerType macro_grid_pointer(macroGridName);
  // refine the grid 'starting_refinement_level' times:
  macro_grid_pointer->globalRefine(coarse_grid_level_);
  //! ---- tools ----
  L2Error< MsFEMTraits::DiscreteFunctionType > l2error;

  //! ---------------------------- grid parts ----------------------------------------------
  // grid part for the global function space, required for MsFEM-macro-problem
  MsFEMTraits::GridPartType gridPart(*macro_grid_pointer);
  MsFEMTraits::GridType& grid = gridPart.grid();
  //! --------------------------------------------------------------------------------------

  // coarse grid
  MsFEMTraits::GridPointerType macro_grid_pointer_coarse(macroGridName);
  macro_grid_pointer_coarse->globalRefine(coarse_grid_level_);
  MsFEMTraits::GridPartType gridPart_coarse(*macro_grid_pointer_coarse);
  MsFEMTraits::GridType& grid_coarse = gridPart_coarse.grid();

  grid.globalRefine(total_refinement_level_ - coarse_grid_level_);

  //! ------------------------- discrete function spaces -----------------------------------
  // the global-problem function space:
  MsFEMTraits::DiscreteFunctionSpaceType discreteFunctionSpace(gridPart);
  MsFEMTraits::DiscreteFunctionSpaceType discreteFunctionSpace_coarse(gridPart_coarse);

  //! --------------------------- coefficient functions ------------------------------------

  // defines the matrix A^{\epsilon} in our global problem  - div ( A^{\epsilon}(\nabla u^{\epsilon} ) = f
  const MsFEMTraits::DiffusionType diffusion_op;
  // define (first) source term:
  const MsFEMTraits::FirstSourceType f; // standard source f

  //! ---------------------------- general output parameters ------------------------------
  // general output parameters
  Dune::myDataOutputParameters outputparam;
  data_output(gridPart, discreteFunctionSpace_coarse, outputparam );

  //! ---------------------- solve MsFEM problem ---------------------------
  //! solution vector
  // solution of the standard finite element method
  MsFEMTraits::DiscreteFunctionType msfem_solution("MsFEM Solution", discreteFunctionSpace);
  msfem_solution.clear();

  MsFEMTraits::DiscreteFunctionType coarse_part_msfem_solution("Coarse Part MsFEM Solution", discreteFunctionSpace);
  coarse_part_msfem_solution.clear();

  MsFEMTraits::DiscreteFunctionType fine_part_msfem_solution("Fine Part MsFEM Solution", discreteFunctionSpace);
  fine_part_msfem_solution.clear();

  const int number_of_level_host_entities = grid_coarse.size(0 /*codim*/);

  // number of layers per coarse grid entity T:  U(T) is created by enrichting T with n(T)-layers.
  MsFEMTraits::MacroMicroGridSpecifierType specifier(discreteFunctionSpace_coarse, discreteFunctionSpace);
  for (int i = 0; i < number_of_level_host_entities; i += 1)
  {
    specifier.setLayer(i, number_of_layers_);
  }
  specifier.setOversamplingStrategy( 3 ); //! Important!

  //! create subgrids:
  const bool silence = false;
  {//this scopes subgridlist
    MsFEMTraits::SubGridListType subgrid_list(specifier, silence);

    // just for Dirichlet zero-boundary condition
    Elliptic_Rigorous_MsFEM_Solver< MsFEMTraits::DiscreteFunctionType > msfem_solver(discreteFunctionSpace);
    msfem_solver.solve_dirichlet_zero(diffusion_op, f, specifier, subgrid_list,
                                      coarse_part_msfem_solution, fine_part_msfem_solution, msfem_solution);

    DSC_LOG_INFO << "Solution output for MsFEM Solution." << std::endl;
    solution_output(msfem_solution, coarse_part_msfem_solution,
                    fine_part_msfem_solution, outputparam,
                    total_refinement_level_, coarse_grid_level_);
  }

  //! ---------------------- solve FEM problem with same (fine) resolution ---------------------------
  //! solution vector
  // solution of the standard finite element method
  MsFEMTraits::DiscreteFunctionType fem_solution("FEM Solution", discreteFunctionSpace);
  fem_solution.clear();

  if ( DSC_CONFIG_GET("rigorous_msfem.fem_comparison",false) )
  {
    // just for Dirichlet zero-boundary condition
    const Elliptic_FEM_Solver< MsFEMTraits::DiscreteFunctionType > fem_solver(discreteFunctionSpace);
    fem_solver.solve_dirichlet_zero(diffusion_op, f, fem_solution);

    //! ----------------------------------------------------------------------
    DSC_LOG_INFO << "Data output for FEM Solution." << std::endl;
    //! -------------------------- writing data output FEM Solution ----------

    // ------------- VTK data output for FEM solution --------------
    // create and initialize output class
    MsFEMTraits::IOTupleType fem_solution_series(&fem_solution);
    outputparam.set_prefix("/fem_solution");
    MsFEMTraits::DataOutputType fem_dataoutput(gridPart.grid(), fem_solution_series, outputparam);

    // write data
    fem_dataoutput.writeData( 1.0 /*dummy*/, "fem_solution" );
    // -------------------------------------------------------------
  }

  DSC_LOG_INFO << std::endl << "The L2 errors:" << std::endl << std::endl;
  //! ----------------- compute L2- and H1- errors -------------------
  if (Problem::ModelProblemData::has_exact_solution)
  {

    H1Error< MsFEMTraits::DiscreteFunctionType > h1error;

    const MsFEMTraits::ExactSolutionType u;
    int order_quadrature_rule = 13;

    MsFEMTraits::RangeType msfem_error = l2error.norm< MsFEMTraits::ExactSolutionType >(u,
                                                              msfem_solution,
                                                              order_quadrature_rule );
    DSC_LOG_INFO << "|| u_msfem - u_exact ||_L2 =  " << msfem_error << std::endl << std::endl;

    MsFEMTraits::RangeType h1_msfem_error(0.0);
    h1_msfem_error = h1error.semi_norm< MsFEMTraits::ExactSolutionType >(u, msfem_solution, order_quadrature_rule);
    h1_msfem_error += msfem_error;
    DSC_LOG_INFO << "|| u_msfem - u_exact ||_H1 =  " << h1_msfem_error << std::endl << std::endl;

    if ( DSC_CONFIG_GET("rigorous_msfem.fem_comparison",false) )
    {

      MsFEMTraits::RangeType approx_msfem_error = l2error.norm2< 2* MsFEMTraits::DiscreteFunctionSpaceType::polynomialOrder + 2 >(fem_solution,
                                                                                                      msfem_solution);
      DSC_LOG_INFO << "|| u_msfem - u_fem ||_L2 =  " << approx_msfem_error << std::endl << std::endl;
      H1Norm< MsFEMTraits::GridPartType > h1norm(gridPart);
      MsFEMTraits::RangeType h1_approx_msfem_error = h1norm.distance(fem_solution, msfem_solution);

      DSC_LOG_INFO << "|| u_msfem - u_fem ||_H1 =  " << h1_approx_msfem_error << std::endl << std::endl;


      MsFEMTraits::RangeType fem_error = l2error.norm< MsFEMTraits::ExactSolutionType >(u,
                                                            fem_solution,
                                                            order_quadrature_rule);

      DSC_LOG_INFO << "|| u_fem - u_exact ||_L2 =  " << fem_error << std::endl << std::endl;

      MsFEMTraits::RangeType h1_fem_error(0.0);

      h1_fem_error = h1error.semi_norm< MsFEMTraits::ExactSolutionType >(u, fem_solution, order_quadrature_rule);
      h1_fem_error += fem_error;
      DSC_LOG_INFO << "|| u_fem - u_exact ||_H1 =  " << h1_fem_error << std::endl << std::endl;
    }
  } else if ( DSC_CONFIG_GET("rigorous_msfem.fem_comparison",false) )
  {
    DSC_LOG_ERROR << "Exact solution not available. Errors between MsFEM and FEM approximations for the same fine grid resolution."
                  << std::endl << std::endl;
    MsFEMTraits::RangeType approx_msfem_error = l2error.norm2< 2* MsFEMTraits::DiscreteFunctionSpaceType::polynomialOrder + 2 >(fem_solution,
                                                                                                      msfem_solution);
    DSC_LOG_INFO << "|| u_msfem - u_fem ||_L2 =  " << approx_msfem_error << std::endl << std::endl;
    H1Norm< MsFEMTraits::GridPartType > h1norm(gridPart);
    MsFEMTraits::RangeType h1_approx_msfem_error = h1norm.distance(fem_solution, msfem_solution);

    DSC_LOG_INFO << "|| u_msfem - u_fem ||_H1 =  " << h1_approx_msfem_error << std::endl << std::endl;
  }

} // function algorithm

} //namespace MsFEM {
} //namespace Multiscale {
} //namespace Dune {

#endif // DUNE_MULTISCALE_RIGOROUS_HH
