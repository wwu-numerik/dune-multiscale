#include "common.hh"

// #define ADAPTIVE

#ifndef ADAPTIVE
 #define UNIFORM
#endif


#if HAVE_GRAPE
 #include <dune/grid/io/visual/grapedatadisplay.hh>
#endif

#define PGF

#include <dune/fem/gridpart/common/gridpart.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/space/lagrangespace.hh>
#include <dune/fem/function/adaptivefunction.hh>

//! ----- typedefs for the macro grid and the corresponding discrete space -----
typedef Dune::GridSelector::GridType
GridType;
// Dune::InteriorBorder_Partition or Dune::All_Partition >?
// see:
// http://www.dune-project.org/doc/doxygen/dune-grid-html/group___g_i_related_types.html#ga5b9e8102d7f70f3f4178182629d98b6
typedef Dune::AdaptiveLeafGridPart< GridType /*,Dune::All_Partition*/ > GridPartType;

typedef Dune::GridPtr< GridType > GridPointerType;

typedef Dune::FunctionSpace< double, double, WORLDDIM, 1 > FunctionSpaceType;

// to display data with ParaView:
#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/fem/io/file/dataoutput.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/io/file/datawriter.hh>

// dune-fem includes:
#include <dune/fem/space/common/adaptmanager.hh>
#include <dune/fem/misc/l2error.hh>
#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/misc/h1norm.hh>

//! local (dune-multiscale) includes
#include <dune/multiscale/problems/elliptic_problems/selector.hh>
#include <dune/multiscale/tools/solver/FEM/fem_solver.hh>
#include <dune/multiscale/tools/solver/MsFEM/msfem_localproblems/subgrid-list.hh>
#include <dune/multiscale/tools/solver/MsFEM/msfem_solver.hh>
#include <dune/multiscale/tools/misc/h1error.hh>
#include <dune/multiscale/tools/disc_func_writer/discretefunctionwriter.hh>
#include <dune/multiscale/tools/meanvalue.hh>
#include <dune/multiscale/tools/improved_l2error.hh>
#include <dune/multiscale/tools/errorestimation/MsFEM/elliptic_error_estimator.hh>

//!-----------------------------------------------------------------------------

//! --------- typedefs for the coefficient and data functions ------------------
using namespace Dune;
// type of first source term (right hand side of differential equation or type of 'f')
typedef Problem::FirstSource< FunctionSpaceType > FirstSourceType;

// type of (possibly non-linear) diffusion term (i.e. 'A^{\epsilon}')
typedef Problem::Diffusion< FunctionSpaceType > DiffusionType;

// default type for any missing coefficient function (e.g. advection,...)
typedef Problem::DefaultDummyFunction< FunctionSpaceType > DefaultDummyFunctionType;

// type of exact solution (in general unknown)
typedef Problem::ExactSolution< FunctionSpaceType > ExactSolutionType;
typedef GridFunctionAdapter< ExactSolutionType, GridPartType >
DiscreteExactSolutionType;     // for data output with paraview or grape


//! ----  typedefs for the standard discrete function space (macroscopic) -----
typedef FunctionSpaceType::DomainType DomainType;
//! define the type of elements of the codomain v(\Omega) (typically a subset of \R)
typedef FunctionSpaceType::RangeType RangeType;
//! defines the function space to which the numerical solution belongs to
//! see dune/fem/lagrangebase.hh
typedef LagrangeDiscreteFunctionSpace< FunctionSpaceType, GridPartType, 1 >  // 1 = POLORDER
DiscreteFunctionSpaceType;

typedef AdaptiveDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;
//!-----------------------------------------------------------------------------

//!------------------------- for adaptive grid refinement ---------------------------------
//! For adaption:
//! type of restrict-prolong operator
typedef RestrictProlongDefault< DiscreteFunctionType >
RestrictProlongOperatorType;
//! type of the adaption manager
typedef AdaptationManager< GridType, RestrictProlongOperatorType >
AdaptationManagerType;
//!---------------------------------------------------------------------------------------

typedef MacroMicroGridSpecifier< DiscreteFunctionSpaceType >                          MacroMicroGridSpecifierType;
typedef SubGrid< GridType::dimension, GridType >                                      SubGridType;
typedef SubGridList< DiscreteFunctionType, SubGridType, MacroMicroGridSpecifierType > SubGridListType;

//! -------------------------- MsFEM error estimator ----------------------------
typedef MsFEMErrorEstimator< DiscreteFunctionType,
                             DiffusionType,
                             FirstSourceType,
                             MacroMicroGridSpecifierType,
                             SubGridListType >
MsFEMErrorEstimatorType;
//! -----------------------------------------------------------------------------

#include "msfem_globals.hh"

//! ------------------ typedefs and classes for data output ---------------------
typedef tuple< DiscreteFunctionType* >      IOTupleType;
typedef DataOutput< GridType, IOTupleType > DataOutputType;
// just for the discretized exact solution (in case it is available)
typedef tuple< const DiscreteExactSolutionType* > ExSolIOTupleType;
// just for the discretized exact solution (in case it is available)
typedef DataOutput< GridType, ExSolIOTupleType > ExSolDataOutputType;

#include <dune/multiscale/tools/misc/outputparameter.hh>

//!---------------------------------------------------------------------------------------

//! algorithm
void algorithm(const std::string& macroGridName) {
  DSC_LOG_INFO << "loading dgf: " << macroGridName << std::endl;
  // we might use further grid parameters (depending on the grid type, e.g. Alberta), here we switch to default values
  // for the parameters:
  // create a grid pointer for the DGF file belongig to the macro grid:
  GridPointerType macro_grid_pointer(macroGridName);
  // refine the grid 'starting_refinement_level' times:
  macro_grid_pointer->globalRefine(coarse_grid_level_);
  //! ---- tools ----
  // model problem data
  Problem::ModelProblemData problem_info;
  L2Error< DiscreteFunctionType > l2error;
  H1Error< DiscreteFunctionType > DUNE_UNUSED(h1error);
  // expensive hack to deal with discrete functions, defined on different grids
  Dune::ImprovedL2Error< DiscreteFunctionType > DUNE_UNUSED(impL2error);

  //! ---------------------------- grid parts ----------------------------------------------
  // grid part for the global function space, required for MsFEM-macro-problem
  GridPartType gridPart(*macro_grid_pointer);
  GridType& grid = gridPart.grid();
  //! --------------------------------------------------------------------------------------

  // coarse grid
  GridPointerType macro_grid_pointer_coarse(macroGridName);
  macro_grid_pointer_coarse->globalRefine(coarse_grid_level_);
  GridPartType gridPart_coarse(*macro_grid_pointer_coarse);
  GridType& grid_coarse = gridPart_coarse.grid();

  // strategy for adaptivity:
  #ifdef ADAPTIVE
  typedef GridType::LeafGridView         GridView;
  typedef GridView::Codim< 0 >::Iterator ElementLeafIterator;
  typedef GridType::Traits::LeafIndexSet GridLeafIndexSet;

  if (local_indicators_available_ == true)
  {
    const bool coarse_scale_error_dominant = false;
    const bool fine_scale_error_dominant = false;   // wird noch nicht benoetigt, dass wir diese Verfeinerung uniform regeln
    const bool oversampling_error_dominant = false;   // wird noch nicht benoetigt, dass wir diese Adadption uniform regeln
    // identify the dominant contribution:
    const double average_est_error = total_estimated_H1_error_[loop_number_ - 1] / 6.0;     // 6 contributions

    if ( (total_approximation_error_[loop_number_ - 1] >= average_est_error)
         || (total_fine_grid_jumps_[loop_number_ - 1] >= average_est_error) )
    {
      fine_scale_error_dominant = true;   /* increase fine level resolution by 2 levels */
      total_refinement_level_ += 2;   // 'the fine grid level'
      DSC_LOG_INFO << "Fine scale error identified as being dominant. Decrease the number of global refinements by 2."
                << std::endl;
      DSC_LOG_INFO << "NEW: Refinement Level for (uniform) Fine Grid = " << total_refinement_level_ << std::endl;
      DSC_LOG_INFO << "Note that this means: the fine grid is " << total_refinement_level_ - coarse_grid_level_
                << " refinement levels finer than the coarse grid." << std::endl;
    }

    if ( (total_projection_error_[loop_number_ - 1] >= average_est_error)
         || (total_conservative_flux_jumps_[loop_number_ - 1] >= average_est_error) )
    {
      oversampling_error_dominant = true;   /* increase number of layers by 1 */
      number_of_layers_ += 1;
      DSC_LOG_INFO
      << "Oversampling error identified as being dominant. Increase the number of layers for each subgrid by 5."
      << std::endl;
      DSC_LOG_INFO << "NEW: Number of layers = " << number_of_layers_ << std::endl;
    }

    if ( (total_coarse_residual_[loop_number_ - 1] >= average_est_error)
         || (total_coarse_grid_jumps_[loop_number_ - 1] >= average_est_error) )
    {
      DSC_LOG_INFO << "Coarse scale error identified as being dominant. Start adaptive coarse grid refinment."
                << std::endl;
      coarse_scale_error_dominant = true;   /* mark elementwise for 2 refinments */
    }

    std::vector< RangeType > average_coarse_error_indicator(loop_number_);        // arithmetic average
    for (int l = 0; l < loop_number_; ++l)
    {
      assert( loc_coarse_residual_[l].size() == loc_coarse_grid_jumps_[l].size() );
      average_coarse_error_indicator[l] = 0.0;
      for (size_t i = 0; i < loc_coarse_grid_jumps_[l].size(); ++i)
      {
        average_coarse_error_indicator[l] += loc_coarse_residual_[l][i] + loc_coarse_grid_jumps_[l][i];
      }
      average_coarse_error_indicator[l] = average_coarse_error_indicator[l] / loc_coarse_residual_[l].size();
    }

    // allgemeineren Algorithmus vorgestellt, aber noch nicht implementiert

    if (coarse_scale_error_dominant == true)
    {
      int number_of_refinements = 4;
      // allowed varianve from average ( in percent )
      double variance = 1.1;   // = 110 % of the average
      for (int l = 0; l < loop_number_; ++l)
      {
        const GridLeafIndexSet& gridLeafIndexSet = grid.leafIndexSet();
        GridView gridView = grid.leafView();

        int total_number_of_entities = 0;
        int number_of_marked_entities = 0;
        for (ElementLeafIterator it = gridView.begin< 0 >();
             it != gridView.end< 0 >(); ++it)
        {
          RangeType loc_coarse_error_indicator = loc_coarse_grid_jumps_[l][gridLeafIndexSet.index(*it)]
                                                 + loc_coarse_residual_[l][gridLeafIndexSet.index(*it)];
          total_number_of_entities += 1;

          if (loc_coarse_error_indicator >= variance * average_coarse_error_indicator[l])
          {
            grid.mark(number_of_refinements, *it);
            number_of_marked_entities += 1;
          }
        }

        if ( l == (loop_number_ - 1) )
        {
          DSC_LOG_INFO << number_of_marked_entities << " of " << total_number_of_entities
                    << " coarse grid entities marked for mesh refinement." << std::endl;
        }

        grid.preAdapt();
        grid.adapt();
        grid.postAdapt();

        // coarse grid
        GridView gridView_coarse = grid_coarse.leafView();
        const GridLeafIndexSet& gridLeafIndexSet_coarse = grid_coarse.leafIndexSet();

        for (ElementLeafIterator it = gridView_coarse.begin< 0 >();
             it != gridView_coarse.end< 0 >(); ++it)
        {
          RangeType loc_coarse_error_indicator = loc_coarse_grid_jumps_[l][gridLeafIndexSet_coarse.index(*it)]
                                                 + loc_coarse_residual_[l][gridLeafIndexSet_coarse.index(*it)];

          if (loc_coarse_error_indicator >= variance * average_coarse_error_indicator[l])
          {
            grid_coarse.mark(number_of_refinements, *it);
          }
        }
        grid_coarse.preAdapt();
        grid_coarse.adapt();
        grid_coarse.postAdapt();
      }
    }
  }

  #endif // #ifdef ADAPTIVE

  grid.globalRefine(total_refinement_level_ - coarse_grid_level_);

  //! ------------------------- discrete function spaces -----------------------------------

  // the global-problem function space:
  DiscreteFunctionSpaceType discreteFunctionSpace(gridPart);

  DiscreteFunctionSpaceType discreteFunctionSpace_coarse(gridPart_coarse);

  //! --------------------------------------------------------------------------------------

  //! --------------------------- coefficient functions ------------------------------------

  // defines the matrix A^{\epsilon} in our global problem  - div ( A^{\epsilon}(\nabla u^{\epsilon} ) = f
  DiffusionType diffusion_op;
  // define (first) source term:
  FirstSourceType f; // standard source f


  //! ---------------------------- general output parameters ------------------------------
  // general output parameters
  myDataOutputParameters outputparam;
  outputparam.set_path(path_);
  // sequence stamp
  std::stringstream outstring;
  //! --------------------------------------------------------------------------------------

  //! -------------------------- writing data output Exact Solution ------------------------
  if (Problem::ModelProblemData::has_exact_solution)
  {
    const ExactSolutionType u;
    const DiscreteExactSolutionType discrete_exact_solution("discrete exact solution ", u, gridPart);

    // create and initialize output class
    ExSolIOTupleType exact_solution_series(&discrete_exact_solution);
    outputparam.set_prefix("exact_solution");
    ExSolDataOutputType exactsol_dataoutput(gridPart.grid(), exact_solution_series, outputparam);

    // write data
    outstring << "exact-solution";
    exactsol_dataoutput.writeData( 1.0 /*dummy*/, outstring.str() );
    // clear the std::stringstream:
    outstring.str( std::string() );
    // -------------------------------------------------------
  }
  //! --------------------------------------------------------------------------------------

  //! --------------- writing data output for the coarse grid visualization ------------------
  DiscreteFunctionType coarse_grid_visualization(filename_ + " Visualization of the coarse grid",
                                                 discreteFunctionSpace_coarse);
  coarse_grid_visualization.clear();
  // -------------------------- data output -------------------------
  // create and initialize output class
  IOTupleType coarse_grid_series(&coarse_grid_visualization);
  char coarse_grid_fname[50];
  sprintf(coarse_grid_fname, "coarse_grid_visualization_%d_", loop_number_);
  const std::string coarse_grid_fname_s(coarse_grid_fname);
  outputparam.set_prefix(coarse_grid_fname_s);
  DataOutputType coarse_grid_dataoutput(gridPart_coarse.grid(), coarse_grid_series, outputparam);
  // write data
  outstring << coarse_grid_fname_s;
  coarse_grid_dataoutput.writeData( 1.0 /*dummy*/, outstring.str() );
  // clear the std::stringstream:
  outstring.str( std::string() );
  // -------------------------------------------------------
  //! --------------------------------------------------------------------------------------

  //! ---------------------- solve MsFEM problem ---------------------------
  //! solution vector
  // solution of the standard finite element method
  DiscreteFunctionType msfem_solution(filename_ + " MsFEM Solution", discreteFunctionSpace);
  msfem_solution.clear();

  DiscreteFunctionType coarse_part_msfem_solution(filename_ + " Coarse Part MsFEM Solution", discreteFunctionSpace);
  coarse_part_msfem_solution.clear();

  DiscreteFunctionType fine_part_msfem_solution(filename_ + " Fine Part MsFEM Solution", discreteFunctionSpace);
  fine_part_msfem_solution.clear();

  const int number_of_level_host_entities = grid_coarse.size(0 /*codim*/);
  const int DUNE_UNUSED(coarse_level_fine_level_difference) = grid.maxLevel() - grid_coarse.maxLevel();

  // number of layers per coarse grid entity T:  U(T) is created by enrichting T with n(T)-layers.
  MacroMicroGridSpecifierType specifier(discreteFunctionSpace_coarse, discreteFunctionSpace);
  for (int i = 0; i < number_of_level_host_entities; i += 1)
  {
    specifier.setLayer(i, number_of_layers_);
  }

  //! create subgrids:
  const bool silence = false;
  SubGridListType subgrid_list(specifier, silence);

  // just for Dirichlet zero-boundary condition
  Elliptic_MsFEM_Solver< DiscreteFunctionType > msfem_solver(discreteFunctionSpace, path_);
  msfem_solver.solve_dirichlet_zero(diffusion_op, f, specifier, subgrid_list,
                                    coarse_part_msfem_solution, fine_part_msfem_solution, msfem_solution);

  DSC_LOG_INFO << "Data output for MsFEM Solution." << std::endl;
  //! ----------------- writing data output MsFEM Solution -----------------
  // --------- VTK data output for MsFEM solution --------------------------
  // create and initialize output class
  IOTupleType msfem_solution_series(&msfem_solution);
  #ifdef ADAPTIVE
  char msfem_fname[50];
  sprintf(msfem_fname, "msfem_solution_%d_", loop_number_);
  std::string msfem_fname_s(msfem_fname);
  outputparam.set_prefix(msfem_fname_s);
  DataOutputType msfem_dataoutput(gridPart.grid(), msfem_solution_series, outputparam);
  // write data
  outstring << msfem_fname_s;
  #else // ifdef ADAPTIVE
  outputparam.set_prefix("msfem_solution");
  DataOutputType msfem_dataoutput(gridPart.grid(), msfem_solution_series, outputparam);
  // write data
  outstring << "msfem_solution";
  #endif // ifdef ADAPTIVE

  msfem_dataoutput.writeData( 1.0 /*dummy*/, outstring.str() );
  // clear the std::stringstream:
  outstring.str( std::string() );

  // create and initialize output class
  IOTupleType coarse_msfem_solution_series(&coarse_part_msfem_solution);

  #ifdef ADAPTIVE
  char coarse_msfem_fname[50];
  sprintf(coarse_msfem_fname, "coarse_part_msfem_solution_%d_", loop_number_);
  const std::string coarse_msfem_fname_s(coarse_msfem_fname);
  outputparam.set_prefix(coarse_msfem_fname_s);
  DataOutputType coarse_msfem_dataoutput(gridPart.grid(), coarse_msfem_solution_series, outputparam);
  // write data
  outstring << coarse_msfem_fname_s;
  #else // ifdef ADAPTIVE
  outputparam.set_prefix("coarse_part_msfem_solution");
  DataOutputType coarse_msfem_dataoutput(gridPart.grid(), coarse_msfem_solution_series, outputparam);
  // write data
  outstring << "coarse_part_msfem_solution";
  #endif // ifdef ADAPTIVE

  coarse_msfem_dataoutput.writeData( 1.0 /*dummy*/, outstring.str() );
  // clear the std::stringstream:
  outstring.str( std::string() );

  // create and initialize output class
  IOTupleType fine_msfem_solution_series(&fine_part_msfem_solution);

  #ifdef ADAPTIVE
  char fine_msfem_fname[50];
  sprintf(fine_msfem_fname, "fine_part_msfem_solution_%d_", loop_number_);
  std::string fine_msfem_fname_s(fine_msfem_fname);
  outputparam.set_prefix(fine_msfem_fname_s);
  DataOutputType fine_msfem_dataoutput(gridPart.grid(), fine_msfem_solution_series, outputparam);
  // write data
  outstring << fine_msfem_fname_s;
  #else // ifdef ADAPTIVE
  outputparam.set_prefix("fine_part_msfem_solution");
  DataOutputType fine_msfem_dataoutput(gridPart.grid(), fine_msfem_solution_series, outputparam);
  // write data
  outstring << "fine_msfem_solution";
  #endif // ifdef ADAPTIVE

  fine_msfem_dataoutput.writeData( 1.0 /*dummy*/, outstring.str() );
  // clear the std::stringstream:
  outstring.str( std::string() );

  // ----------------------------------------------------------------------
  // ---------------------- write discrete msfem solution to file ---------
  char fname[50];
  sprintf(fname, "/msfem_solution_discFunc_refLevel_%d_%d", total_refinement_level_, coarse_grid_level_);
  std::string fname_s(fname);

  std::string location = path_ + fname_s;
  DiscreteFunctionWriter dfw( (location).c_str() );
  if (dfw.is_open())
    dfw.append(msfem_solution);
  //! --------------------------------------------------------------------

  //! ------------------------------------------------------------

  DSC_LOG_INFO << std::endl << "The L2 errors:" << std::endl << std::endl;
  //! ----------------- compute L2- and H1- errors -------------------
  if (Problem::ModelProblemData::has_exact_solution)
  {
    const ExactSolutionType u;
    RangeType msfem_error = l2error.norm< ExactSolutionType >(u,
                                                              msfem_solution,
                                                              2 * DiscreteFunctionSpaceType::polynomialOrder + 2);
    DSC_LOG_INFO << "|| u_msfem - u_exact ||_L2 =  " << msfem_error << std::endl << std::endl;

    RangeType h1_msfem_error(0.0);
    h1_msfem_error = h1error.semi_norm< ExactSolutionType >(u, msfem_solution);
    h1_msfem_error += msfem_error;
    DSC_LOG_INFO << "|| u_msfem - u_exact ||_H1 =  " << h1_msfem_error << std::endl << std::endl;
  }

  // error estimation
  RangeType total_estimated_H1_error(0.0);

  MsFEMErrorEstimatorType estimator(discreteFunctionSpace, specifier, subgrid_list, diffusion_op, f, path_);
  total_estimated_H1_error = estimator.adaptive_refinement(grid_coarse,
                                                           msfem_solution,
                                                           coarse_part_msfem_solution,
                                                           fine_part_msfem_solution);
  {//intentional scope
    std::vector<RangeVectorVector*> locals = {{ &loc_coarse_residual_, &loc_coarse_grid_jumps_, &loc_projection_error_, &loc_conservative_flux_jumps_, &loc_approximation_error_, &loc_fine_grid_jumps_}};
    std::vector<RangeVector*> totals = {{&total_coarse_residual_, &total_projection_error_, &total_coarse_grid_jumps_, &total_conservative_flux_jumps_, &total_approximation_error_, &total_fine_grid_jumps_ }};
    assert(locals.size() == totals.size());
    for(auto loc : locals) (*loc)[loop_number_] = RangeVector(specifier.getNumOfCoarseEntities(),0.0);
    for(auto total : totals) (*total)[loop_number_] = 0.0;

    for (int m = 0; m < specifier.getNumOfCoarseEntities(); ++m) {
      loc_coarse_residual_[loop_number_][m] = specifier.get_loc_coarse_residual(m);
      loc_coarse_grid_jumps_[loop_number_][m] = specifier.get_loc_coarse_grid_jumps(m);
      loc_projection_error_[loop_number_][m] = specifier.get_loc_projection_error(m);
      loc_conservative_flux_jumps_[loop_number_][m] = specifier.get_loc_conservative_flux_jumps(m);
      loc_approximation_error_[loop_number_][m] = specifier.get_loc_approximation_error(m);
      loc_fine_grid_jumps_[loop_number_][m] = specifier.get_loc_fine_grid_jumps(m);

      for (size_t i = 0; i < totals.size(); ++i) (*totals[i])[loop_number_] += std::pow((*locals[i])[loop_number_][m], 2.0);
    }
    total_estimated_H1_error_[loop_number_] = 0.0;
    for(auto total : totals) {
      (*total)[loop_number_] = std::sqrt((*total)[loop_number_]);
      total_estimated_H1_error_[loop_number_] += (*total)[loop_number_];
    }
    local_indicators_available_ = true;
  }

  #ifdef ADAPTIVE
  if (total_estimated_H1_error <= error_tolerance_)
    repeat_algorithm_ = false;
  #else // ifdef ADAPTIVE
  repeat_algorithm_ = false;
  #endif // ifdef ADAPTIVE

  //! ---------------------- solve FEM problem ---------------------------
  //! solution vector
  // solution of the standard finite element method
  DiscreteFunctionType fem_solution(filename_ + " FEM Solution", discreteFunctionSpace);
  fem_solution.clear();

  // just for Dirichlet zero-boundary condition
  const Elliptic_FEM_Solver< DiscreteFunctionType > fem_solver(discreteFunctionSpace);
  fem_solver.solve_dirichlet_zero(diffusion_op, f, fem_solution);

  //! ----------------------------------------------------------------------
  DSC_LOG_INFO << "Data output for FEM Solution." << std::endl;
  //! -------------------------- writing data output FEM Solution ----------

  // ------------- VTK data output for FEM solution --------------
  // create and initialize output class
  IOTupleType fem_solution_series(&fem_solution);
  outputparam.set_prefix("fem_solution");
  DataOutputType fem_dataoutput(gridPart.grid(), fem_solution_series, outputparam);

  // write data
  outstring << "fem_solution";
  fem_dataoutput.writeData( 1.0 /*dummy*/, outstring.str() );
  // clear the std::stringstream:
  outstring.str( std::string() );

  // -------------------------------------------------------------

  // ------------- write discrete fem solution to file -----------
  char fem_fname[50];
  sprintf(fem_fname, "/fem_solution_discFunc_refLevel_%d", total_refinement_level_);
  std::string fem_fname_s(fem_fname);

  std::string fine_location = path_ + fem_fname_s;
  DiscreteFunctionWriter fem_dfw( (fine_location).c_str() );
  if (fem_dfw.is_open())
    fem_dfw.append(fem_solution);

  DSC_LOG_INFO << std::endl << "The L2 errors:" << std::endl << std::endl;
  //! ----------------- compute L2- and H1- errors -------------------
  if (Problem::ModelProblemData::has_exact_solution)
  {
    const ExactSolutionType u;
    RangeType fem_error = l2error.norm< ExactSolutionType >(u,
                                                            fem_solution,
                                                            2 * DiscreteFunctionSpaceType::polynomialOrder + 2);

    DSC_LOG_INFO << "|| u_fem - u_exact ||_L2 =  " << fem_error << std::endl << std::endl;

    RangeType h1_fem_error(0.0);
    h1_fem_error = h1error.semi_norm< ExactSolutionType >(u, fem_solution);
    h1_fem_error += fem_error;
    DSC_LOG_INFO << "|| u_fem - u_exact ||_H1 =  " << h1_fem_error << std::endl << std::endl;
  } else {
    DSC_LOG_ERROR << "Exact solution not available. Use fine-scale FEM-approximation as a reference solution."
                  << std::endl << std::endl;
    RangeType approx_msfem_error = l2error.norm2< 2* DiscreteFunctionSpaceType::polynomialOrder + 2 >(fem_solution,
                                                                                                      msfem_solution);
    DSC_LOG_INFO << "|| u_msfem - u_fem ||_L2 =  " << approx_msfem_error << std::endl << std::endl;
    H1Norm< GridPartType > h1norm(gridPart);
    RangeType h1_approx_msfem_error = h1norm.distance(fem_solution, msfem_solution);

    DSC_LOG_INFO << "|| u_msfem - u_fem ||_H1 =  " << h1_approx_msfem_error << std::endl << std::endl;
  }
  //! -------------------------------------------------------

} // function algorithm

int main(int argc, char** argv) {
  try {
    init(argc, argv);
    namespace DSC = Dune::Stuff::Common;
    //!TODO include base in config
    path_ = std::string("data/MsFEM/") + DSC::Parameter::Config().get("global.datadir", "data");

    // generate directories for data output
    DSC::Filesystem::testCreateDirectory(path_);

    std::string save_filename = path_ + "/problem-info.txt";
    DSC_LOG_INFO << "Data will be saved under: " << save_filename << std::endl;

    const int start_level = DSC::Parameter::Config().get("grid.start_level", 4);
    coarse_grid_level_ = DSC::Parameter::Config().get( "grid.coarse_level", 4, DSC::Parameter::ValidateLess< int >(start_level) );
    number_of_layers_ = DSC::Parameter::Config().get("grid.oversampling_layers", 4);

    #ifdef ADAPTIVE
    error_tolerance_ = DSC::Parameter::Config().get("problem.error_tolerance", 1e-6);
    #endif   // ifdef ADAPTIVE

    // data for the model problem; the information manager
    // (see 'problem_specification.hh' for details)
    const Problem::ModelProblemData info(filename_);

    // refinement_level denotes the (starting) grid refinement level for the global problem, i.e. it describes 'H'
    total_refinement_level_
      = DSC::Parameter::Config().get( "grid.total_refinement", 4, DSC::Parameter::ValidateLess< int >(coarse_grid_level_) );

    // name of the grid file that describes the macro-grid:
    const std::string macroGridName = info.getMacroGridFile();

    DSC_LOG_INFO << "Error File for Elliptic Model Problem " << info.get_Number_of_Model_Problem()
              << " with epsilon = " << info.getEpsilon() << "." << std::endl << std::endl;
    #ifdef UNIFORM
    DSC_LOG_INFO << "Use MsFEM with an uniform computation, i.e.:" << std::endl;
    DSC_LOG_INFO << "Uniformly refined coarse and fine mesh and" << std::endl;
    DSC_LOG_INFO << "the same number of layers for each (oversampled) local grid computation." << std::endl << std::endl;
    DSC_LOG_INFO << "Computations were made for:" << std::endl << std::endl;
    DSC_LOG_INFO << "Refinement Level for (uniform) Fine Grid = " << total_refinement_level_ << std::endl;
    DSC_LOG_INFO << "Refinement Level for (uniform) Coarse Grid = " << coarse_grid_level_ << std::endl;
    DSC_LOG_INFO << "Number of layers for oversampling = " << number_of_layers_ << std::endl;
    DSC_LOG_INFO << std::endl << std::endl;
    #else   // ifdef UNIFORM
    DSC_LOG_INFO << "Use MsFEM with an adaptive computation, i.e.:" << std::endl;
    DSC_LOG_INFO << "Starting with a uniformly refined coarse and fine mesh and" << std::endl;
    DSC_LOG_INFO << "the same number of layers for each (oversampled) local grid computation." << std::endl << std::endl;
    DSC_LOG_INFO << "Error tolerance = " << error_tolerance_ << std::endl << std::endl;
    DSC_LOG_INFO << "Computations were made for:" << std::endl << std::endl;
    DSC_LOG_INFO << "(Starting) Refinement Level for (uniform) Fine Grid = " << total_refinement_level_ << std::endl;
    DSC_LOG_INFO << "(Starting) Refinement Level for (uniform) Coarse Grid = " << coarse_grid_level_ << std::endl;
    DSC_LOG_INFO << "(Starting) Number of layers for oversampling = " << number_of_layers_ << std::endl;
    DSC_LOG_INFO << std::endl << std::endl;
    #endif   // ifdef UNIFORM

    loop_number_ = 0;
    while (repeat_algorithm_ == true)
    {
      #ifdef ADAPTIVE
      DSC_LOG_INFO << "------------------ run " << loop_number_ + 1 << " --------------------" << std::endl << std::endl;
      #endif
      algorithm(macroGridName);
      #ifdef ADAPTIVE
      DSC_LOG_INFO << std::endl << std::endl;
      DSC_LOG_INFO << "---------------------------------------------" << std::endl;
      loop_number_ += 1;
      #endif   // ifdef ADAPTIVE
    }
    // the reference problem generaly has a 'refinement_difference_for_referenceproblem' higher resolution than the
    //normal
    // macro problem

    const long double cpu_time = clock() / double(CLOCKS_PER_SEC);
    DSC_LOG_INFO << "Total runtime of the program: " << cpu_time << "s" << std::endl;
    return 0;
  } catch (Dune::Exception& e) {
    std::cerr << e.what() << std::endl;
  }
  return 1;
} // main
