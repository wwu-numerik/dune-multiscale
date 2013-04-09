// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "common.hh"

#include <dune/multiscale/msfem/algorithm.hh>

int main(int argc, char** argv) {
  try {
    init(argc, argv);
    using namespace Dune::Multiscale::MsFEM;

    //!TODO include base in config
    DSC_PROFILER.startTiming("msfem.all");

    const std::string datadir = DSC_CONFIG_GET("global.datadir", "data/");

    // generate directories for data output
    DSC::testCreateDirectory(datadir);

    DSC_LOG_INFO << boost::format("Data will be saved under: %s\nLogs will be saved under: %s/%s/ms.log.log\n")
                            % datadir % datadir % DSC_CONFIG_GET("logging.dir", "log");

    // syntax: info_from_par_file / default  / validation of the value

    // coarse_grid_level denotes the (starting) grid refinement level for the global coarse scale problem, i.e. it describes 'H'
    int coarse_grid_level_ = DSC_CONFIG_GETV( "msfem.coarse_grid_level", 4, DSC::ValidateLess< int >( -1 ) );

    // syntax: info_from_par_file / default
    int number_of_layers_ = DSC_CONFIG_GET("msfem.oversampling_layers", 4);

    switch ( DSC_CONFIG_GET( "msfem.oversampling_strategy", 1 ) )
    {
      case 1: break;
      case 2: break;
      default: DUNE_THROW(Dune::InvalidStateException, "Oversampling Strategy must be 1 or 2.");
    }
      //if (!( (DSC_CONFIG_GET( "msfem.oversampling_strategy", 1 ) == 1) || (DSC_CONFIG_GET( "msfem.oversampling_strategy", 1 ) == 2) ))
     
    
    // data for the model problem; the information manager
    // (see 'problem_specification.hh' for details)
    const Problem::ModelProblemData info;

    // total_refinement_level denotes the (starting) grid refinement level for the global fine scale problem, i.e. it describes 'h'
    int total_refinement_level_
      = DSC_CONFIG_GETV( "msfem.fine_grid_level", 4, DSC::ValidateLess< int >(coarse_grid_level_-1) );

    // name of the grid file that describes the macro-grid:
    const std::string macroGridName = info.getMacroGridFile();

    DSC_LOG_INFO << "Error File for Elliptic Model Problem " << Dune::Stuff::Common::getTypename(info)
              << " with epsilon = " << DSC_CONFIG_GET("problem.epsilon", 1.0f) << "." << std::endl << std::endl;
    if (DSC_CONFIG_GET("msfem.uniform", true)) {
      if ( DSC_CONFIG_GET("msfem.petrov_galerkin", true) )
        DSC_LOG_INFO << "Use MsFEM in Petrov-Galerkin formulation with an uniform computation, i.e.:" << std::endl;
      else
        DSC_LOG_INFO << "Use MsFEM in classical (symmetric) formulation with an uniform computation, i.e.:" << std::endl;      
      DSC_LOG_INFO << "Uniformly refined coarse and fine mesh and" << std::endl;
      DSC_LOG_INFO << "the same number of layers for each (oversampled) local grid computation." << std::endl << std::endl;
      DSC_LOG_INFO << "Computations were made for:" << std::endl << std::endl;
      DSC_LOG_INFO << "Refinement Level for (uniform) Fine Grid = " << total_refinement_level_ << std::endl;
      DSC_LOG_INFO << "Refinement Level for (uniform) Coarse Grid = " << coarse_grid_level_ << std::endl;
      DSC_LOG_INFO << "Oversampling Strategy = " << DSC_CONFIG_GET( "msfem.oversampling_strategy", 1 ) << std::endl;
      DSC_LOG_INFO << "Number of layers for oversampling = " << number_of_layers_ << std::endl;
      if ( DSC_CONFIG_GET("msfem.fem_comparison",false) )
       { DSC_LOG_INFO << std::endl << "Comparison with standard FEM computation on the MsFEM Fine Grid, i.e. on Refinement Level " << total_refinement_level_ << std::endl; }
      DSC_LOG_INFO << std::endl << std::endl;
    } else {
      if ( DSC_CONFIG_GET("msfem.petrov_galerkin", true) )
        DSC_LOG_INFO << "Use MsFEM in Petrov-Galerkin formulation with an adaptive computation, i.e.:" << std::endl;
      else
        DSC_LOG_INFO << "Use MsFEM in classical (symmetric) formulation with an adaptive computation, i.e.:" << std::endl;  
      DSC_LOG_INFO << "Starting with a uniformly refined coarse and fine mesh and" << std::endl;
      DSC_LOG_INFO << "the same number of layers for each (oversampled) local grid computation." << std::endl << std::endl;
      DSC_LOG_INFO << "Error tolerance = " << DSC_CONFIG_GET("msfem.error_tolerance", 1e-6) << std::endl << std::endl;
      DSC_LOG_INFO << "Computations were made for:" << std::endl << std::endl;
      DSC_LOG_INFO << "(Starting) Refinement Level for (uniform) Fine Grid = " << total_refinement_level_ << std::endl;
      DSC_LOG_INFO << "(Starting) Refinement Level for (uniform) Coarse Grid = " << coarse_grid_level_ << std::endl;
      DSC_LOG_INFO << "Oversampling Strategy = " << DSC_CONFIG_GET( "msfem.oversampling_strategy", 1 ) << std::endl;
      DSC_LOG_INFO << "(Starting) Number of layers for oversampling = " << number_of_layers_ << std::endl;
      if ( DSC_CONFIG_GET("msfem.fem_comparison",false) )
      { DSC_LOG_INFO << std::endl << "Comparison with a standard FEM computation on the MsFEM Fine Grid." << std::endl; }
      DSC_LOG_INFO << std::endl << std::endl;
    }

    //! ---------------------- local error indicators --------------------------------
    // ----- local error indicators (for each coarse grid element T) -------------
    const int max_loop_number = 10;
    // local coarse residual, i.e. H ||f||_{L^2(T)}
    MsFEMTraits::RangeVectorVector loc_coarse_residual_(max_loop_number);
    // local coarse grid jumps (contribute to the total coarse residual)
    MsFEMTraits::RangeVectorVector loc_coarse_grid_jumps_(max_loop_number);
    // local projection error (we project to get a globaly continous approximation)
    MsFEMTraits::RangeVectorVector loc_projection_error_(max_loop_number);
    // local jump in the conservative flux
    MsFEMTraits::RangeVectorVector loc_conservative_flux_jumps_(max_loop_number);
    // local approximation error
    MsFEMTraits::RangeVectorVector loc_approximation_error_(max_loop_number);
    // local sum over the fine grid jumps (for a fixed subgrid that cooresponds with a coarse entity T)
    MsFEMTraits::RangeVectorVector loc_fine_grid_jumps_(max_loop_number);

    MsFEMTraits::RangeVector total_coarse_residual_(max_loop_number);
    MsFEMTraits::RangeVector total_projection_error_(max_loop_number);
    MsFEMTraits::RangeVector total_coarse_grid_jumps_(max_loop_number);
    MsFEMTraits::RangeVector total_conservative_flux_jumps_(max_loop_number);
    MsFEMTraits::RangeVector total_approximation_error_(max_loop_number);
    MsFEMTraits::RangeVector total_fine_grid_jumps_(max_loop_number);
    MsFEMTraits::RangeVector total_estimated_H1_error_(max_loop_number);

    //! TODO put these into something like a named tuple/class
    std::vector<MsFEMTraits::RangeVectorVector*> locals = {{ &loc_coarse_residual_, &loc_coarse_grid_jumps_,
                                                             &loc_projection_error_, &loc_conservative_flux_jumps_,
                                                             &loc_approximation_error_, &loc_fine_grid_jumps_}};
    std::vector<MsFEMTraits::RangeVector*> totals = {{&total_coarse_residual_, &total_projection_error_,
                                                      &total_coarse_grid_jumps_, &total_conservative_flux_jumps_,
                                                      &total_approximation_error_, &total_fine_grid_jumps_ }};

    unsigned int loop_number = 0;
    while (algorithm(macroGridName, loop_number++, total_refinement_level_, coarse_grid_level_,
                     number_of_layers_, locals, totals, total_estimated_H1_error_))
    {}

    // the reference problem generaly has a 'refinement_difference_for_referenceproblem' higher resolution than the
    //normal
    // macro problem

    const auto cpu_time = DSC_PROFILER.stopTiming("msfem.all");
    DSC_LOG_INFO << "Total runtime of the program: " << cpu_time << "ms" << std::endl;
    DSC_PROFILER.outputTimings("profiler");
    return 0;
  } catch (Dune::Exception& e) {
    std::cerr << e.what() << std::endl;
  }
  return 1;
} // main
