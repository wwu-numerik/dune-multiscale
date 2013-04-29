// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "common.hh"

#include <dune/multiscale/common/traits.hh>
#include <dune/multiscale/hmm/algorithm.hh>

// the solutions of the cell problems are always determined in a pre-process

#include <dune/multiscale/tools/errorestimation/HMM/elliptic_error_estimator.hh>

//! (very restrictive) homogenizer
#include <dune/multiscale/tools/homogenizer/elliptic_analytical_homogenizer.hh>
#include <dune/multiscale/tools/homogenizer/elliptic_homogenizer.hh>

// all the multiscale code requires an access to the 'ModelProblemData' class
// (typically defined in the specification file for the model problem)

//! local (dune-multiscale) includes
#include <dune/multiscale/tools/righthandside_assembler.hh>
#include <dune/multiscale/tools/disc_func_writer/discretefunctionwriter.hh>
#include <dune/multiscale/tools/meanvalue.hh>

void check_config();

int main(int argc, char** argv) {
  try {
    init(argc, argv);
    check_config();

    using namespace Dune::Multiscale;
    using namespace Dune::Multiscale::HMM;

    DSC_PROFILER.startTiming("total_cpu");

    const std::string path = DSC_CONFIG_GET("global.datadir", "data/");

    // generate directories for data output
    DSC::testCreateDirectory(path);

    if ( DSC_CONFIG_GET("problem.stochastic_pertubation", false)) {
      //! Do we want to force the algorithm to come to an end?
      // (was auch immer der Grund war, dass das Programm zuvor endlos lange weiter gelaufen ist. z.B. Tolerenzen nicht
      // erreicht etc.)
      DSC_CONFIG.set("problem.force_end", true);
    }

    // name of the error file in which the data will be saved
    std::string filename_;
    const Problem::ModelProblemData info;

    // man koennte hier noch den genauen Iterationsschritt in den Namen mit einfliessen lassen:
    // (vorlauefig sollte diese Variante aber reichen)
    const std::string save_filename = DSC_CONFIG_GET("RESUME_TO_BROKEN_COMPUTATION", false)
                                      ? std::string(path + "problem-info-resumed-computation.txt")
                                      : std::string(path + "/logdata/ms.log.log");
    DSC_LOG_INFO << "LOG FILE " << std::endl << std::endl;
    DSC_LOG_INFO << "Data will be saved under: " << save_filename << std::endl;

    // refinement_level denotes the (starting) grid refinement level for the global problem, i.e. it describes 'H'
    const int refinement_level_macrogrid_ = DSC_CONFIG_GET("hmm.coarse_grid_level", 4);
    // grid refinement level for solving the cell problems, i.e. it describes 'h':
    const int refinement_level_cellgrid = DSC_CONFIG_GET("hmm.cell_grid_level", 4);
    // (macro) grid refinement level for which the reference problem was solved
    int refinement_level_referenceprob_ = 0;
    if (DSC_CONFIG_GET("problem.reference_solution", false)) // if there is a refernce solution
      refinement_level_referenceprob_ = DSC_CONFIG_GET("problem.rs_grid_level", 0);
    // (typically this is a very high level so that we get a very fine grid)

    // name of the grid file that describes the macro-grid:
    const std::string macroGridName = info.getMacroGridFile();
    DSC_LOG_INFO << "loading dgf: " << macroGridName << std::endl;

    // we might use further grid parameters (depending on the grid type, e.g. Alberta), here we switch to default values
    // for the parameters:

    // create a grid pointer for the DGF file belongig to the macro grid:
    CommonTraits::GridPointerType macro_grid_pointer(macroGridName);
    // refine the grid 'starting_refinement_level' times:
    macro_grid_pointer->globalRefine(refinement_level_macrogrid_);

    // create a finer GridPart for either the homogenized or the fine-scale problem.
    // this shall be used to compute an approximation of the exact solution.
    CommonTraits::GridPointerType fine_macro_grid_pointer(macroGridName);
    // refine the grid 'starting_refinement_level_reference' times:
    fine_macro_grid_pointer->globalRefine(refinement_level_referenceprob_);

    // after transformation, the cell problems are problems on the 0-centered unit cube [-½,½]²:
    // THIS GRID IS FIXED!
    const boost::filesystem::path unitCubeName = "../dune/multiscale/grids/cell_grids/unit_cube_0_centered.dgf";
    // --> the 0-centered unit cube, i.e. [-1/2,1/2]^2
    // note that the centering is fundamentaly important for the implementation. Do NOT change it to e.g. [0,1]^2!!!
    // to solve the cell problems, we always need a periodic gridPart.
    // Here it is always the unit cube that needs to be used (after transformation, cell problems are always formulated
    // on such a grid )
    Dune::GridPtr< CommonTraits::GridType > periodic_grid_pointer(unitCubeName.string());
    periodic_grid_pointer->globalRefine(refinement_level_cellgrid);

    algorithm(macro_grid_pointer, fine_macro_grid_pointer,
              periodic_grid_pointer, filename_);
    // the reference problem generaly has a 'refinement_difference_for_referenceproblem' higher resolution than the
    //normal
    // macro problem

    const auto cpu_time = DSC_PROFILER.stopTiming("total_cpu") / 1000.f;
    DSC_LOG_INFO << "Total runtime of the program: " << cpu_time << "ms" << std::endl;
    return 0;
  } catch (Dune::Exception& e) {
    std::cerr << e.what() << std::endl;
  }
  return 1;
} // main

//! \brief sanity checks for our configuration
void check_config()
{
  if ( !DSC_CONFIG_GET("hmm.error_estimation", false) && DSC_CONFIG_GET("hmm.adaptivity", false) )
    DUNE_THROW(Dune::InvalidStateException, "Error estimation must be activated to use adaptivity.");

  if (DSC_CONFIG_GET("RESUME_TO_BROKEN_COMPUTATION", false)) {
    DSC_CONFIG.set("HMM_NEWTON_ITERATION_STEP", 2);
  }
  else {
    DSC_CONFIG.set("HMM_NEWTON_ITERATION_STEP", 0);
  }

}
