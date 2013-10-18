#include <config.h>
// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

// Implementation of the Local Orthogonal Decomposition Method (LODM)

#include <dune/multiscale/common/main_init.hh>
#include <dune/multiscale/lod/algorithm.hh>
#include <dune/multiscale/problems/selector.hh>

int main(int argc, char** argv) {
  try {
    using namespace Dune::Multiscale;
    using namespace Dune::Multiscale::MsFEM;
    init(argc, argv);

    //!TODO include base in config
    DSC_PROFILER.startTiming("msfem.all");

    const std::string datadir = DSC_CONFIG_GET("global.datadir", "data/");

    // generate directories for data output
    DSC::testCreateDirectory(datadir);

    DSC_LOG_INFO << boost::format("Data will be saved under: %s\nLogs will be saved under: %s/%s/ms.log.log\n") %
                        datadir % datadir % DSC_CONFIG_GET("logging.dir", "log");

    // syntax: info_from_par_file / default  / validation of the value

    // coarse_grid_level denotes the (starting) grid refinement level for the global coarse scale problem, i.e. it
    // describes 'H'
    int coarse_grid_level_ = DSC_CONFIG_GETV("rigorous_msfem.coarse_grid_level", 4, DSC::ValidateLess<int>(-1));

    // syntax: info_from_par_file / default
    int number_of_layers_ = DSC_CONFIG_GET("rigorous_msfem.oversampling_layers", 4);

    if (!((DSC_CONFIG_GET("rigorous_msfem.oversampling_strategy", "Clement") == "Clement") ||
          (DSC_CONFIG_GET("rigorous_msfem.oversampling_strategy", "Clement") == "Lagrange"))) {
      DUNE_THROW(Dune::InvalidStateException, "Oversampling Strategy must be 'Lagrange' or 'Clement'.");
    }

    // total_refinement_level denotes the (starting) grid refinement level for the global fine scale problem, i.e. it
    // describes 'h'
    int total_refinement_level_ =
        DSC_CONFIG_GETV("rigorous_msfem.fine_grid_level", 4, DSC::ValidateLess<int>(coarse_grid_level_ - 1));

    // name of the grid file that describes the macro-grid:
    auto info = Problem::getModelData();
    const std::string macroGridName = info->getMacroGridFile();

    DSC_LOG_INFO << "Error File for Elliptic Model Problem " << DSC::getTypename(*info)
                 << " with epsilon = " << DSC_CONFIG_GET("problem.epsilon", 1.0f) << "." << std::endl << std::endl;
    if (DSC_CONFIG_GET("lod.petrov_galerkin", false))
      DSC_LOG_INFO << "Use Local Orthogonal Decomposition (LOD) Method in Petrov-Galerkin formulation with an uniform "
                      "computation, i.e.:" << std::endl;
    else
      DSC_LOG_INFO << "Use Local Orthogonal Decomposition (LOD) Method in classical (symmetric) formulation with an "
                      "uniform computation, i.e.:" << std::endl;
    DSC_LOG_INFO << "Uniformly refined coarse and fine mesh and" << std::endl;
    DSC_LOG_INFO << "the same number of layers for each (oversampled) local grid computation." << std::endl
                 << std::endl;
    DSC_LOG_INFO << "Computations were made for:" << std::endl << std::endl;
    DSC_LOG_INFO << "Refinement Level for (uniform) Fine Grid = " << total_refinement_level_ << std::endl;
    DSC_LOG_INFO << "Refinement Level for (uniform) Coarse Grid = " << coarse_grid_level_ << std::endl;
    // DSC_LOG_INFO << "Oversampling Strategy = " << DSC_CONFIG_GET( "rigorous_msfem.oversampling_strategy", 1 ) <<
    // std::endl;
    DSC_LOG_INFO << "Number of layers for oversampling = " << number_of_layers_ << std::endl;
    DSC_LOG_INFO << "Oversampling Strategy = " << DSC_CONFIG_GET("rigorous_msfem.oversampling_strategy", "Clement")
                 << std::endl;
    if (DSC_CONFIG_GET("rigorous_msfem.fem_comparison", false)) {
      DSC_LOG_INFO << std::endl
                   << "Comparison with standard FEM computation on the MsFEM Fine Grid, i.e. on Refinement Level "
                   << total_refinement_level_ << std::endl;
    }
    DSC_LOG_INFO << std::endl << std::endl;

    algorithm(macroGridName, total_refinement_level_, coarse_grid_level_, number_of_layers_);

    // the reference problem generaly has a 'refinement_difference_for_referenceproblem' higher resolution than the
    // normal
    // macro problem

    const auto cpu_time = DSC_PROFILER.stopTiming("msfem.all");
    DSC_LOG_INFO << "Total runtime of the program: " << cpu_time << "ms" << std::endl;
    DSC_PROFILER.outputTimings("profiler");
    return 0;
  }
  catch (Dune::Exception& e) {
    std::cerr << e.what() << std::endl;
  }
  return 1;
} // main
