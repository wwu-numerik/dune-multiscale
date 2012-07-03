#include "common.hh"
#include "hmm_config.hh"

#include <dune/multiscale/hmm/algorithm.hh>

// we only use error estimation, if the solutions of the cell problems have been determined in a pre-process. Otherwise
// it is far too expensive!
#ifndef AD_HOC_COMPUTATION
 #include <dune/multiscale/tools/errorestimation/HMM/elliptic_error_estimator.hh>
#endif

// ! (very restrictive) homogenizer
#ifdef LINEAR_PROBLEM
 #include <dune/multiscale/tools/homogenizer/elliptic_analytical_homogenizer.hh>
 #include <dune/multiscale/tools/homogenizer/elliptic_homogenizer.hh>
#else // ifdef LINEAR_PROBLEM
      // dummy (does not work, since identical to HMM assembler)
// was deleted #include <dune/multiscale/tools/homogenizer/nonlinear_elliptic_homogenizer.hh>
#endif // ifdef LINEAR_PROBLEM

// ! NOTE: All the multiscale code requires an access to the 'ModelProblemData' class (typically defined in
// problem_specification.hh), which provides us with information about epsilon, delta, etc.
// ! HMM Assembler, Error Estimator, ... they all hark back to 'ModelProblemData'. Probably there is a better solution,
// but for me, it works perfectly.
namespace Multiscale {
// parameters for the current realization of the HMM
class HMMParameters
{
public:
  // Constructor for ModelProblemData
  inline explicit HMMParameters()
  {}

  // do you want to save the solutions of the cell problems for a later usage?
  // (e.g. for error estimation and adaptivity)
  inline bool save_cell_problems() const {
    return true;
  }
};
}

// ! local (dune-multiscale) includes
#include <dune/multiscale/tools/assembler/righthandside_assembler.hh>
#include <dune/multiscale/tools/disc_func_writer/discretefunctionwriter.hh>
#include <dune/multiscale/tools/solver/HMM/reconstruction_manager/elliptic/reconstructionintegrator.hh>
#include <dune/multiscale/tools/meanvalue.hh>

int main(int argc, char** argv) {
  try {
    init(argc, argv);

    std::string path = std::string("data/HMM/") + Stuff::Config().get("global.datadir", "data");
    // generate directories for data output
    Stuff::testCreateDirectory(path);

    // ! -------------------------------- important variables ---------------------------------

    // name of the error file in which the data will be saved
    std::string filename_;

    double epsilon_; // 'epsilon' in for instance A^{epsilon}(t,x) = A(t,x/epsilon)
    double epsilon_est_; // estimated epsilon in case epsilon is unknown
    double delta_; // edge length of the cells in the cell proplems,

    int refinement_level_macrogrid_;
    int refinement_level_referenceprob_;

    #ifdef RESUME_TO_BROKEN_COMPUTATION
    // man koennte hier noch den genauen Iterationsschritt in den Namen mit einfliessen lassen:
    // (vorlauefig sollte diese Variante aber reichen)
    std::string save_filename = path + "/problem-info-resumed-computation.txt";
    #else // ifdef RESUME_TO_BROKEN_COMPUTATION
    std::string save_filename = path + "/problem-info.txt";
    #endif // ifdef RESUME_TO_BROKEN_COMPUTATION
    std::cout << "Data will be saved under: " << save_filename << std::endl;
    // data for the model problem; the information manager
    // (see 'problem_specification.hh' for details)
    Problem::ModelProblemData info(path);
    // epsilon is specified in ModelProblemData, which is specified in problem_specification.hh
    epsilon_ = info.getEpsilon();
    // estimated epsilon (specified in ModelProblemData)
    epsilon_est_ = info.getEpsilonEstimated();
    // edge length of the cells in the cells, belonging to the cell problems
    // note that (delta/epsilon_est) needs to be a positive integer!
    delta_ = info.getDelta();

    // refinement_level denotes the (starting) grid refinement level for the global problem, i.e. it describes 'H'
    refinement_level_macrogrid_ = Stuff::Config().get("grid.refinement_level_macrogrid", 0);

    // grid refinement level for solving the cell problems, i.e. it describes 'h':
    const int refinement_level_cellgrid = Stuff::Config().get("grid.refinement_level_cellgrid", 1);

    #ifdef ADAPTIVE
    error_tolerance_ = Stuff::Config().get("problem.error_tolerance", 1e-6);
    #endif // ifdef ADAPTIVE

    // (starting) grid refinement level for solving the reference problem
    refinement_level_referenceprob_ = info.getRefinementLevelReferenceProblem();
    // in general: for the homogenized case = 11 and for the high resolution case = 14
    // Note that this depends on the model problem!
    #ifndef FINE_SCALE_REFERENCE
    refinement_level_referenceprob_ = 8;
    #endif

    // how many times finer do we solve the reference problem (it is either a homogenized problem or the exact problem
    // with a fine-scale resolution)
    int refinement_difference_for_referenceproblem = refinement_level_referenceprob_ - refinement_level_macrogrid_;
    // name of the grid file that describes the macro-grid:
    std::string macroGridName;
    info.getMacroGridFile(macroGridName);
    std::cout << "loading dgf: " << macroGridName << std::endl;

    // we might use further grid parameters (depending on the grid type, e.g. Alberta), here we switch to default values
    // for the parameters:

    // create a grid pointer for the DGF file belongig to the macro grid:
    HMM::GridPointerType macro_grid_pointer(macroGridName);
    // refine the grid 'starting_refinement_level' times:
    macro_grid_pointer->globalRefine(refinement_level_macrogrid_);

    // create a finer GridPart for either the homogenized or the fine-scale problem.
    // this shall be used to compute an approximation of the exact solution.
    HMM::GridPointerType fine_macro_grid_pointer(macroGridName);
    // refine the grid 'starting_refinement_level_reference' times:
    fine_macro_grid_pointer->globalRefine(refinement_level_referenceprob_);

    // after transformation, the cell problems are problems on the 0-centered unit cube [-½,½]²:
    std::string UnitCubeName("../dune/multiscale/grids/cell_grids/unit_cube_0_centered.dgf");     // --> the 0-centered
                                                                                                  // unit cube, i.e.
                                                                                                  // [-1/2,1/2]^2
    // note that the centering is fundamentaly important for the implementation. Do NOT change it to e.g. [0,1]^2!!!
    // to solve the cell problems, we always need a periodic gridPart.
    // Here it is always the unit cube that needs to be used (after transformation, cell problems are always formulated
    //on
    // such a grid )
    Dune::GridPtr< HMM::GridType > periodic_grid_pointer(UnitCubeName);
    periodic_grid_pointer->globalRefine(refinement_level_cellgrid);

    // to save all information in a file
    std::ofstream data_file( (save_filename).c_str() );
    if ( data_file.is_open() )
    {
      data_file << "Error File for Elliptic Model Problem " << info.get_Number_of_Model_Problem() << "." << std::endl
                << std::endl;
      #ifdef LINEAR_PROBLEM
      data_file << "Problem is declared as being LINEAR." << std::endl;
      #else
      data_file << "Problem is declared as being NONLINEAR." << std::endl;
      #endif // ifdef LINEAR_PROBLEM
      #ifdef EXACTSOLUTION_AVAILABLE
      data_file << "Exact solution is available." << std::endl << std::endl;
      #else
      data_file << "Exact solution is not available." << std::endl << std::endl;
      #endif // ifdef EXACTSOLUTION_AVAILABLE
      data_file << "Computations were made for:" << std::endl << std::endl;
      data_file << "Refinement Level for (uniform) Macro Grid = " << refinement_level_macrogrid_ << std::endl;
      data_file << "Refinement Level for Periodic Micro Grid = " << refinement_level_cellgrid << std::endl << std::endl;
      #ifdef TFR
      data_file << "We use TFR-HMM (HMM with test function reconstruction)." << std::endl;
      #else
      data_file << "We use HMM without test function reconstruction (NO TFR)." << std::endl;
      #endif // ifdef TFR
      #ifdef AD_HOC_COMPUTATION
      data_file << "Cell problems are solved ad hoc (where required)." << std::endl << std::endl;
      #else
      data_file << "Cell problems are solved and saved (in a pre-process)." << std::endl << std::endl;
      #ifdef ERRORESTIMATION
      data_file << "Error estimation activated!" << std::endl << std::endl;
      #endif
      #endif // ifdef AD_HOC_COMPUTATION
      data_file << "Epsilon = " << epsilon_ << std::endl;
      data_file << "Estimated Epsilon = " << epsilon_est_ << std::endl;
      data_file << "Delta (edge length of cell-cube) = " << delta_ << std::endl;
      #ifdef STOCHASTIC_PERTURBATION
      data_file << std::endl << "Stochastic perturbation added. Variance = " << VARIANCE << std::endl;
      #endif
      #ifdef ADAPTIVE
      data_file << std::endl << "Adaptive computation. Global error tolerance for program abort = "
                << error_tolerance_ << std::endl;
      #endif // ifdef ADAPTIVE
      data_file << std::endl << std::endl;
    }

    algorithm(UnitCubeName,
              macro_grid_pointer,
              fine_macro_grid_pointer,
              periodic_grid_pointer,
              refinement_difference_for_referenceproblem,
              data_file,
              filename_);
    // the reference problem generaly has a 'refinement_difference_for_referenceproblem' higher resolution than the
    //normal
    // macro problem

    long double cpu_time = clock();
    cpu_time = cpu_time / CLOCKS_PER_SEC;
    std::cout << "Total runtime of the program: " << cpu_time << "s" << std::endl;

    if ( data_file.is_open() )
    {
      data_file << "Total runtime of the program: " << cpu_time << "s" << std::endl;
    }
    data_file.close();
    return 0;
  } catch (Dune::Exception& e) {
    std::cerr << e.what() << std::endl;
  }
  return 1;
} // main
