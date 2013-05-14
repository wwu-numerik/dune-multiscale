// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "algorithm.hh"

#include <dune/fem/gridpart/common/gridpart.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/space/lagrangespace.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/grid/common/gridinfo.hh>

// to display data with ParaView:
#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/fem/io/file/dataoutput.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/io/file/datawriter.hh>
#include <dune/fem/space/common/adaptmanager.hh>
#include <dune/fem/misc/l2error.hh>
#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/misc/h1norm.hh>

#include <dune/stuff/common/filesystem.hh>
//! local (dune-multiscale) includes
#include <dune/multiscale/problems/elliptic_problems/selector.hh>
#include <dune/multiscale/msfem/fem_solver.hh>
#include <dune/multiscale/msfem/localproblems/subgrid-list.hh>
#include <dune/multiscale/msfem/msfem_solver.hh>
#include <dune/multiscale/msfem/fem_solver.hh>
#include <dune/multiscale/tools/misc/h1error.hh>
#include <dune/multiscale/hmm/cell_problem_numbering.hh>
#include <dune/multiscale/tools/meanvalue.hh>
#include <dune/multiscale/tools/improved_l2error.hh>
#include <dune/multiscale/msfem/msfem_elliptic_error_estimator.hh>
#include <dune/multiscale/tools/misc/outputparameter.hh>
#include <dune/multiscale/msfem/msfem_grid_specifier.hh>
#include <dune/multiscale/problems/elliptic_problems/selector.hh>

#include <dune/multiscale/msfem/msfem_traits.hh>
#include <dune/multiscale/common/traits.hh>

namespace Dune {
namespace Multiscale {
namespace MsFEM {

//! \TODO docme
void adapt(CommonTraits::GridType& grid,
           CommonTraits::GridType& grid_coarse,
           const int loop_number,
           int& total_refinement_level_,
           int& coarse_grid_level_,
           int& number_of_layers_,
           const std::vector<CommonTraits::RangeVectorVector*>& locals,
           const std::vector<CommonTraits::RangeVector*>& totals,
           const CommonTraits::RangeVector& total_estimated_H1_error_)
{
  using namespace Dune;
  typedef CommonTraits::GridType::LeafGridView         GridView;
  typedef GridView::Codim< 0 >::Iterator ElementLeafIterator;
  typedef CommonTraits::GridType::Traits::LeafIndexSet GridLeafIndexSet;

  bool coarse_scale_error_dominant = false;
  // identify the dominant contribution:
  const double average_est_error = total_estimated_H1_error_[loop_number - 1] / 6.0;     // 6 contributions

  const auto& total_approximation_error_ = *totals[4];
  const auto& total_fine_grid_jumps_ = *totals[5];
  if ( (total_approximation_error_[loop_number - 1] >= average_est_error)
       || (total_fine_grid_jumps_[loop_number - 1] >= average_est_error) )
  {
    total_refinement_level_ += 2;   // 'the fine grid level'
    DSC_LOG_INFO << "Fine scale error identified as being dominant. Decrease the number of global refinements by 2."
              << std::endl;
    DSC_LOG_INFO << "NEW: Refinement Level for (uniform) Fine Grid = " << total_refinement_level_ << std::endl;
    DSC_LOG_INFO << "Note that this means: the fine grid is " << total_refinement_level_ - coarse_grid_level_
              << " refinement levels finer than the coarse grid." << std::endl;
  }

  const auto& total_projection_error_ = *totals[1];
  const auto& total_conservative_flux_jumps_ = *totals[3];
  if ( (total_projection_error_[loop_number - 1] >= average_est_error)
       || (total_conservative_flux_jumps_[loop_number - 1] >= average_est_error) )
  {
    number_of_layers_ += 1;
    DSC_LOG_INFO
    << "Oversampling error identified as being dominant. Increase the number of layers for each subgrid by 5."
    << std::endl;
    DSC_LOG_INFO << "NEW: Number of layers = " << number_of_layers_ << std::endl;
  }

  const auto& total_coarse_residual_ = *totals[0];
  const auto& total_coarse_grid_jumps_ = *totals[2];
  if ( (total_coarse_residual_[loop_number - 1] >= average_est_error)
       || (total_coarse_grid_jumps_[loop_number - 1] >= average_est_error) )
  {
    DSC_LOG_INFO << "Coarse scale error identified as being dominant. Start adaptive coarse grid refinment."
              << std::endl;
    coarse_scale_error_dominant = true;   /* mark elementwise for 2 refinments */
  }

  const auto& loc_coarse_residual_ = *locals[0];
  const auto& loc_coarse_grid_jumps_ = *totals[1];
  std::vector< CommonTraits::RangeType > average_coarse_error_indicator(loop_number);        // arithmetic average
  for (int l = 0; l < loop_number; ++l)
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

  if (coarse_scale_error_dominant)
  {
    int number_of_refinements = 4;
    // allowed varianve from average ( in percent )
    double variance = 1.1;   // = 110 % of the average
    for (int l = 0; l < loop_number; ++l)
    {
      const auto& gridLeafIndexSet = grid.leafIndexSet();
      auto gridView = grid.leafView();

      int total_number_of_entities = 0;
      int number_of_marked_entities = 0;
      for (ElementLeafIterator it = gridView.begin< 0 >();
           it != gridView.end< 0 >(); ++it)
      {
        CommonTraits::RangeType loc_coarse_error_indicator = loc_coarse_grid_jumps_[l][gridLeafIndexSet.index(*it)]
                                               + loc_coarse_residual_[l][gridLeafIndexSet.index(*it)];
        total_number_of_entities += 1;

        if (loc_coarse_error_indicator >= variance * average_coarse_error_indicator[l])
        {
          grid.mark(number_of_refinements, *it);
          number_of_marked_entities += 1;
        }
      }

      if ( l == (loop_number - 1) )
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
        CommonTraits::RangeType loc_coarse_error_indicator = loc_coarse_grid_jumps_[l][gridLeafIndexSet_coarse.index(*it)]
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

//! \TODO docme
void solution_output(const CommonTraits::DiscreteFunctionType& msfem_solution,
                     const CommonTraits::DiscreteFunctionType& coarse_part_msfem_solution,
                     const CommonTraits::DiscreteFunctionType& fine_part_msfem_solution,
                     Dune::Multiscale::OutputParameters& outputparam,
                     const int loop_number,
                     int& total_refinement_level_,
                     int& coarse_grid_level_)
{
  using namespace Dune;
  //! ----------------- writing data output MsFEM Solution -----------------
  // --------- VTK data output for MsFEM solution --------------------------
  // create and initialize output class
  CommonTraits::IOTupleType msfem_solution_series(&msfem_solution);
  const auto& gridPart = msfem_solution.space().gridPart();
  std::string outstring;
  if (DSC_CONFIG_GET("adaptive", false)) {
    const std::string msfem_fname_s = (boost::format("/msfem_solution_%d_") % loop_number).str();
    outputparam.set_prefix(msfem_fname_s);
    outstring = msfem_fname_s;
  } else {
    outputparam.set_prefix("/msfem_solution");
    outstring = "msfem_solution";
  }

  CommonTraits::DataOutputType msfem_dataoutput(gridPart.grid(), msfem_solution_series, outputparam);
  msfem_dataoutput.writeData( 1.0 /*dummy*/, outstring );

  // create and initialize output class
  CommonTraits::IOTupleType coarse_msfem_solution_series(&coarse_part_msfem_solution);

  if (DSC_CONFIG_GET("adaptive", false)) {
    const std::string coarse_msfem_fname_s = (boost::format("/coarse_part_msfem_solution_%d_") % loop_number).str();
    outputparam.set_prefix(coarse_msfem_fname_s);
    outstring = coarse_msfem_fname_s;
  } else {
    outputparam.set_prefix("/coarse_part_msfem_solution");
    outstring = "coarse_part_msfem_solution";
  }

  CommonTraits::DataOutputType coarse_msfem_dataoutput(gridPart.grid(), coarse_msfem_solution_series, outputparam);
  coarse_msfem_dataoutput.writeData( 1.0 /*dummy*/, outstring );

  // create and initialize output class
  CommonTraits::IOTupleType fine_msfem_solution_series(&fine_part_msfem_solution);

  if (DSC_CONFIG_GET("adaptive", false)) {
    const std::string fine_msfem_fname_s = (boost::format("/fine_part_msfem_solution_%d_") % loop_number).str();
    outputparam.set_prefix(fine_msfem_fname_s);
    // write data
    outstring = fine_msfem_fname_s;
  } else {
    outputparam.set_prefix("/fine_part_msfem_solution");
    // write data
    outstring = "fine_msfem_solution";
  }

  CommonTraits::DataOutputType fine_msfem_dataoutput(gridPart.grid(), fine_msfem_solution_series, outputparam);
  fine_msfem_dataoutput.writeData( 1.0 /*dummy*/, outstring);

  // ----------------------------------------------------------------------
  // ---------------------- write discrete msfem solution to file ---------
  const std::string location = (boost::format("msfem_solution_discFunc_refLevel_%d_%d")
                                %  total_refinement_level_ % coarse_grid_level_).str();
  DiscreteFunctionWriter(location).append(msfem_solution);
  //! --------------------------------------------------------------------
}

//! \TODO docme
void data_output(const CommonTraits::GridPartType& gridPart,
                 const CommonTraits::DiscreteFunctionSpaceType& discreteFunctionSpace_coarse,
                 Dune::Multiscale::OutputParameters& outputparam,
                 const int loop_number)
{
  using namespace Dune;
  // sequence stamp
  //! --------------------------------------------------------------------------------------

  //! -------------------------- writing data output Exact Solution ------------------------
  if (Problem::ModelProblemData::has_exact_solution)
  {
    const CommonTraits::ExactSolutionType u;
    const CommonTraits::DiscreteExactSolutionType discrete_exact_solution("discrete exact solution ", u, gridPart);
    // create and initialize output class
    CommonTraits::ExSolIOTupleType exact_solution_series(&discrete_exact_solution);
    outputparam.set_prefix("/exact_solution");
    CommonTraits::ExSolDataOutputType exactsol_dataoutput(gridPart.grid(), exact_solution_series, outputparam);
    // write data
    exactsol_dataoutput.writeData( 1.0 /*dummy*/, "exact-solution" );
    // -------------------------------------------------------
  }
  //! --------------------------------------------------------------------------------------

  //! --------------- writing data output for the coarse grid visualization ------------------
  CommonTraits::DiscreteFunctionType coarse_grid_visualization("Visualization of the coarse grid",
                                                 discreteFunctionSpace_coarse);
  coarse_grid_visualization.clear();
  // -------------------------- data output -------------------------
  // create and initialize output class
  CommonTraits::IOTupleType coarse_grid_series(&coarse_grid_visualization);

  const auto coarse_grid_fname = (boost::format("/coarse_grid_visualization_%d_") % loop_number).str();
  outputparam.set_prefix(coarse_grid_fname);
  CommonTraits::DataOutputType coarse_grid_dataoutput(discreteFunctionSpace_coarse.gridPart().grid(), coarse_grid_series, outputparam);
  // write data
  coarse_grid_dataoutput.writeData( 1.0 /*dummy*/, coarse_grid_fname );
  // -------------------------------------------------------
  //! --------------------------------------------------------------------------------------
}

//! \TODO docme
bool error_estimation(const CommonTraits::DiscreteFunctionType& msfem_solution,
                      const CommonTraits::DiscreteFunctionType& coarse_part_msfem_solution,
                      const CommonTraits::DiscreteFunctionType& fine_part_msfem_solution,
                      MsFEMTraits::MsFEMErrorEstimatorType& estimator,
                      MsFEMTraits::MacroMicroGridSpecifierType& specifier,
                      const int loop_number,
                      std::vector<CommonTraits::RangeVectorVector*>& locals,
                      std::vector<CommonTraits::RangeVector*>& totals,
                      CommonTraits::RangeVector& total_estimated_H1_error_)
{
  using namespace Dune;

  CommonTraits::RangeType total_estimated_H1_error(0.0);

  // error estimation
  total_estimated_H1_error = estimator.adaptive_refinement(msfem_solution,
                                                           coarse_part_msfem_solution,
                                                           fine_part_msfem_solution);
  {//intentional scope
    assert(locals.size() == totals.size());
    for(auto loc : locals) (*loc)[loop_number] = CommonTraits::RangeVector(specifier.getNumOfCoarseEntities(),0.0);
    for(auto total : totals) (*total)[loop_number] = 0.0;

    for (int m = 0; m < specifier.getNumOfCoarseEntities(); ++m) {
      (*locals[0])[loop_number][m] = specifier.get_loc_coarse_residual(m);
      (*locals[1])[loop_number][m] = specifier.get_loc_coarse_grid_jumps(m);
      (*locals[2])[loop_number][m] = specifier.get_loc_projection_error(m);
      (*locals[3])[loop_number][m] = specifier.get_loc_conservative_flux_jumps(m);
      (*locals[4])[loop_number][m] = specifier.get_loc_approximation_error(m);
      (*locals[5])[loop_number][m] = specifier.get_loc_fine_grid_jumps(m);

      for (size_t i = 0; i < totals.size(); ++i) (*totals[i])[loop_number] += std::pow((*locals[i])[loop_number][m], 2.0);
    }
    total_estimated_H1_error_[loop_number] = 0.0;
    for(auto total : totals) {
      (*total)[loop_number] = std::sqrt((*total)[loop_number]);
      total_estimated_H1_error_[loop_number] += (*total)[loop_number];
    }
  }

  return DSC_CONFIG_GET("adaptive", false)
      ? total_estimated_H1_error > DSC_CONFIG_GET("msfem.error_tolerance", 1e-6) : false;
}

//! algorithm
bool algorithm(const std::string& macroGridName,
               const int loop_number,
               int& total_refinement_level_,
               int& coarse_grid_level_,
               int& number_of_layers_,
               std::vector<CommonTraits::RangeVectorVector*>& locals,
               std::vector<CommonTraits::RangeVector*>& totals,
               CommonTraits::RangeVector& total_estimated_H1_error_) {
  using namespace Dune;
  bool local_indicators_available_ = false;

  if (DSC_CONFIG_GET("adaptive", false))
    DSC_LOG_INFO << "------------------ run " << loop_number + 1 << " --------------------" << std::endl << std::endl;

  DSC_LOG_INFO << "loading dgf: " << macroGridName << std::endl;
  // we might use further grid parameters (depending on the grid type, e.g. Alberta), here we switch to default values
  // for the parameters:
  // create a grid pointer for the DGF file belongig to the macro grid:
  CommonTraits::GridPointerType macro_grid_pointer(macroGridName);
  // refine the grid 'starting_refinement_level' times:
  macro_grid_pointer->globalRefine(coarse_grid_level_);
  //! ---- tools ----
  L2Error< CommonTraits::DiscreteFunctionType > l2error;

  CommonTraits::GridType& grid = *macro_grid_pointer;
  grid.loadBalance();
  gridinfo(grid);
  //! ---------------------------- grid parts ----------------------------------------------
  // grid part for the global function space, required for MsFEM-macro-problem
  CommonTraits::GridPartType gridPart(grid);

  //! --------------------------------------------------------------------------------------

  // coarse grid
  CommonTraits::GridPointerType macro_grid_pointer_coarse(macroGridName);
  macro_grid_pointer_coarse->globalRefine(coarse_grid_level_);
  CommonTraits::GridPartType gridPart_coarse(*macro_grid_pointer_coarse);
  CommonTraits::GridType& grid_coarse = gridPart_coarse.grid();

  // strategy for adaptivity:
  if (DSC_CONFIG_GET("adaptive", false) && local_indicators_available_)
    adapt(grid, grid_coarse, loop_number, total_refinement_level_, coarse_grid_level_,
          number_of_layers_, locals, totals, total_estimated_H1_error_);

  grid.globalRefine(total_refinement_level_ - coarse_grid_level_);

  //! ------------------------- discrete function spaces -----------------------------------
  // the global-problem function space:
  CommonTraits::DiscreteFunctionSpaceType discreteFunctionSpace(gridPart);
  CommonTraits::DiscreteFunctionSpaceType discreteFunctionSpace_coarse(gridPart_coarse);

  //! --------------------------- coefficient functions ------------------------------------

  // defines the matrix A^{\epsilon} in our global problem  - div ( A^{\epsilon}(\nabla u^{\epsilon} ) = f
  const CommonTraits::DiffusionType diffusion_op;
  // define (first) source term:
  const CommonTraits::FirstSourceType f; // standard source f

  //! ---------------------------- general output parameters ------------------------------
  // general output parameters
  Dune::Multiscale::OutputParameters outputparam;
  data_output(gridPart, discreteFunctionSpace_coarse, outputparam, loop_number);

  //! ---------------------- solve MsFEM problem ---------------------------
  //! solution vector
  // solution of the standard finite element method
  CommonTraits::DiscreteFunctionType msfem_solution("MsFEM Solution", discreteFunctionSpace);
  msfem_solution.clear();

  CommonTraits::DiscreteFunctionType coarse_part_msfem_solution("Coarse Part MsFEM Solution", discreteFunctionSpace);
  coarse_part_msfem_solution.clear();

  CommonTraits::DiscreteFunctionType fine_part_msfem_solution("Fine Part MsFEM Solution", discreteFunctionSpace);
  fine_part_msfem_solution.clear();

  const int number_of_level_host_entities = grid_coarse.size(0 /*codim*/);

  // number of layers per coarse grid entity T:  U(T) is created by enrichting T with n(T)-layers.
  MsFEMTraits::MacroMicroGridSpecifierType specifier(discreteFunctionSpace_coarse, discreteFunctionSpace);
  for (int i = 0; i < number_of_level_host_entities; ++i)
  {
    specifier.setNoOfLayers(i, number_of_layers_);
  }
  specifier.setOversamplingStrategy( DSC_CONFIG_GET( "msfem.oversampling_strategy", 1 ) );

  //default to false, error_estimation might change it
  bool repeat_algorithm = false;

  //! create subgrids:
  const bool silence = false;
  {//this scopes subgridlist
    MsFEMTraits::SubGridListType subgrid_list(specifier, silence);

    // just for Dirichlet zero-boundary condition
    Elliptic_MsFEM_Solver msfem_solver(discreteFunctionSpace);
    msfem_solver.solve_dirichlet_zero(diffusion_op, f, specifier, subgrid_list,
                                      coarse_part_msfem_solution, fine_part_msfem_solution, msfem_solution);

    DSC_LOG_INFO << "Solution output for MsFEM Solution." << std::endl;
    solution_output(msfem_solution, coarse_part_msfem_solution,
                    fine_part_msfem_solution, outputparam, loop_number,
                    total_refinement_level_, coarse_grid_level_);

    // error estimation
    if ( DSC_CONFIG_GET("msfem.error_estimation", 0) ) {
      MsFEMTraits::MsFEMErrorEstimatorType estimator(discreteFunctionSpace, specifier, subgrid_list, diffusion_op, f);
      error_estimation(msfem_solution, coarse_part_msfem_solution,
                       fine_part_msfem_solution, estimator, specifier, loop_number,
                       locals, totals, total_estimated_H1_error_);
      local_indicators_available_ = true;
    }
  }


  //! ---------------------- solve FEM problem with same (fine) resolution ---------------------------
  //! solution vector
  // solution of the standard finite element method
  CommonTraits::DiscreteFunctionType fem_solution("FEM Solution", discreteFunctionSpace);
  fem_solution.clear();

  if ( DSC_CONFIG_GET("msfem.fem_comparison",false) )
  {
    // just for Dirichlet zero-boundary condition
    typedef Dune::Multiscale::Elliptic_FEM_Solver KOP;
    const KOP fem_solver(discreteFunctionSpace);
    fem_solver.solve_dirichlet_zero(diffusion_op, f, fem_solution);

    //! ----------------------------------------------------------------------
    DSC_LOG_INFO << "Data output for FEM Solution." << std::endl;
    //! -------------------------- writing data output FEM Solution ----------

    // ------------- VTK data output for FEM solution --------------
    // create and initialize output class
    CommonTraits::IOTupleType fem_solution_series(&fem_solution);
    outputparam.set_prefix("/fem_solution");
    CommonTraits::DataOutputType fem_dataoutput(gridPart.grid(), fem_solution_series, outputparam);

    // write data
    fem_dataoutput.writeData( 1.0 /*dummy*/, "fem_solution" );
    // -------------------------------------------------------------
  }

  DSC_LOG_INFO << std::endl << "The L2 errors:" << std::endl << std::endl;
  //! ----------------- compute L2- and H1- errors -------------------
  if (Problem::ModelProblemData::has_exact_solution)
  {

    H1Error< CommonTraits::DiscreteFunctionType > h1error;

    const CommonTraits::ExactSolutionType u;
    int order_quadrature_rule = 13;

    CommonTraits::RangeType msfem_error = l2error.norm< CommonTraits::ExactSolutionType >(u,
                                                              msfem_solution,
                                                              order_quadrature_rule );
    DSC_LOG_INFO << "|| u_msfem - u_exact ||_L2 =  " << msfem_error << std::endl << std::endl;

    CommonTraits::RangeType h1_msfem_error(0.0);
    h1_msfem_error = h1error.semi_norm< CommonTraits::ExactSolutionType >(u, msfem_solution, order_quadrature_rule);
    h1_msfem_error += msfem_error;
    DSC_LOG_INFO << "|| u_msfem - u_exact ||_H1 =  " << h1_msfem_error << std::endl << std::endl;

    if ( DSC_CONFIG_GET("msfem.fem_comparison",false) )
    {
      CommonTraits::RangeType fem_error = l2error.norm< CommonTraits::ExactSolutionType >(u,
                                                            fem_solution,
                                                            order_quadrature_rule );

      DSC_LOG_INFO << "|| u_fem - u_exact ||_L2 =  " << fem_error << std::endl << std::endl;

      CommonTraits::RangeType h1_fem_error(0.0);

      h1_fem_error = h1error.semi_norm< CommonTraits::ExactSolutionType >(u, fem_solution, order_quadrature_rule);
      h1_fem_error += fem_error;
      DSC_LOG_INFO << "|| u_fem - u_exact ||_H1 =  " << h1_fem_error << std::endl << std::endl;
    }
  } else if ( DSC_CONFIG_GET("msfem.fem_comparison",false) )
  {
    DSC_LOG_ERROR << "Exact solution not available. Errors between MsFEM and FEM approximations for the same fine grid resolution."
                  << std::endl << std::endl;
    CommonTraits::RangeType approx_msfem_error = l2error.norm2< 2* CommonTraits::DiscreteFunctionSpaceType::polynomialOrder + 2 >(fem_solution,
                                                                                                      msfem_solution);
    DSC_LOG_INFO << "|| u_msfem - u_fem ||_L2 =  " << approx_msfem_error << std::endl << std::endl;
    H1Norm< CommonTraits::GridPartType > h1norm(gridPart);
    CommonTraits::RangeType h1_approx_msfem_error = h1norm.distance(fem_solution, msfem_solution);

    DSC_LOG_INFO << "|| u_msfem - u_fem ||_H1 =  " << h1_approx_msfem_error << std::endl << std::endl;
  }
  //! -------------------------------------------------------
  if (DSC_CONFIG_GET("adaptive", false)) {
    DSC_LOG_INFO << std::endl << std::endl;
    DSC_LOG_INFO << "---------------------------------------------" << std::endl;
  }
  return repeat_algorithm;
} // function algorithm

} //namespace MsFEM {
} //namespace Multiscale {
} //namespace Dune {
