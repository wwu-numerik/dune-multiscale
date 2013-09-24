// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include <dune/multiscale/common/main_init.hh>
#include <dune/multiscale/common/error_container.hh>
#include <dune/multiscale/msfem/algorithm.hh>
#include <dune/multiscale/problems/selector.hh>

// for rusage
#include <sys/resource.h>

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
    DSC_LOG_INFO_0 << boost::format("Data will be saved under: %s\nLogs will be saved under: %s/%s/ms.log.log\n") %
                          datadir % datadir % DSC_CONFIG_GET("logging.dir", "log");

    // syntax: info_from_par_file / default  / validation of the value

    // coarse_grid_level denotes the (starting) grid refinement level for the global coarse scale problem, i.e. it
    // describes 'H'
    int coarse_grid_level_ = DSC_CONFIG_GETV("msfem.coarse_grid_level", 4, DSC::ValidateLess<int>(-1));

    // syntax: info_from_par_file / default
    int number_of_layers_ = DSC_CONFIG_GET("msfem.oversampling_layers", 4);

    switch (DSC_CONFIG_GET("msfem.oversampling_strategy", 1)) {
      case 1:
        break;
      case 2:
        break;
      default:
        DUNE_THROW(Dune::InvalidStateException, "Oversampling Strategy must be 1 or 2.");
    }

    // data for the model problem; the information manager
    // (see 'problem_specification.hh' for details)
    auto info_ptr = Problem::getModelData();
    const auto& info = *info_ptr;

    // total_refinement_level denotes the (starting) grid refinement level for the global fine scale problem, i.e. it
    // describes 'h'
    int total_refinement_level_ =
        DSC_CONFIG_GETV("msfem.fine_grid_level", 4, DSC::ValidateLess<int>(coarse_grid_level_ - 1));

    // name of the grid file that describes the macro-grid:
    const std::string macroGridName = info.getMacroGridFile();

    DSC_LOG_INFO_0 << "Error File for Elliptic Model Problem " << Dune::Stuff::Common::getTypename(info)
                   << " with epsilon = " << DSC_CONFIG_GET("problem.epsilon", 1.0f) << "." << std::endl << std::endl;
    if (DSC_CONFIG_GET("msfem.uniform", true)) {
      if (DSC_CONFIG_GET("msfem.petrov_galerkin", true))
        DSC_LOG_INFO_0 << "Use MsFEM in Petrov-Galerkin formulation with an uniform computation, i.e.:" << std::endl;
      else
        DSC_LOG_INFO_0 << "Use MsFEM in classical (symmetric) formulation with an uniform computation, i.e.:"
                       << std::endl;
      DSC_LOG_INFO_0 << "Uniformly refined coarse and fine mesh and" << std::endl;
      DSC_LOG_INFO_0 << "the same number of layers for each (oversampled) local grid computation." << std::endl
                     << std::endl;
      DSC_LOG_INFO_0 << "Computations were made for:" << std::endl << std::endl;
      DSC_LOG_INFO_0 << "Refinement Level for (uniform) Fine Grid = " << total_refinement_level_ << std::endl;
      DSC_LOG_INFO_0 << "Refinement Level for (uniform) Coarse Grid = " << coarse_grid_level_ << std::endl;
      DSC_LOG_INFO_0 << "Oversampling Strategy = " << DSC_CONFIG_GET("msfem.oversampling_strategy", 1) << std::endl;
      DSC_LOG_INFO_0 << "Number of layers for oversampling = " << number_of_layers_ << std::endl;
      if (DSC_CONFIG_GET("msfem.fem_comparison", false)) {
        DSC_LOG_INFO_0 << std::endl
                       << "Comparison with standard FEM computation on the MsFEM Fine Grid, i.e. on Refinement Level "
                       << total_refinement_level_ << std::endl;
      }
      DSC_LOG_INFO_0 << std::endl << std::endl;
    } else {
      if (DSC_CONFIG_GET("msfem.petrov_galerkin", true))
        DSC_LOG_INFO_0 << "Use MsFEM in Petrov-Galerkin formulation with an adaptive computation, i.e.:" << std::endl;
      else
        DSC_LOG_INFO_0 << "Use MsFEM in classical (symmetric) formulation with an adaptive computation, i.e.:"
                       << std::endl;
      DSC_LOG_INFO_0 << "Starting with a uniformly refined coarse and fine mesh and" << std::endl;
      DSC_LOG_INFO_0 << "the same number of layers for each (oversampled) local grid computation." << std::endl
                     << std::endl;
      DSC_LOG_INFO_0 << "Error tolerance = " << DSC_CONFIG_GET("msfem.error_tolerance", 1e-6) << std::endl << std::endl;
      DSC_LOG_INFO_0 << "Computations were made for:" << std::endl << std::endl;
      DSC_LOG_INFO_0 << "(Starting) Refinement Level for (uniform) Fine Grid = " << total_refinement_level_
                     << std::endl;
      DSC_LOG_INFO_0 << "(Starting) Refinement Level for (uniform) Coarse Grid = " << coarse_grid_level_ << std::endl;
      DSC_LOG_INFO_0 << "Oversampling Strategy = " << DSC_CONFIG_GET("msfem.oversampling_strategy", 1) << std::endl;
      DSC_LOG_INFO_0 << "(Starting) Number of layers for oversampling = " << number_of_layers_ << std::endl;
      if (DSC_CONFIG_GET("msfem.fem_comparison", false)) {
        DSC_LOG_INFO_0 << std::endl << "Comparison with a standard FEM computation on the MsFEM Fine Grid."
                       << std::endl;
      }
      DSC_LOG_INFO_0 << std::endl << std::endl;
    }

    //! ---------------------- local error indicators --------------------------------
    // ----- local error indicators (for each coarse grid element T) -------------
    const int max_loop_number = 10;
    ErrorContainer errors(max_loop_number);

    unsigned int loop_number = 0;
    while (algorithm(macroGridName, loop_number++, total_refinement_level_, coarse_grid_level_, number_of_layers_,
                     errors.locals, errors.totals, errors.total_estimated_H1_error_)) {
    }

    // the reference problem generaly has a 'refinement_difference_for_referenceproblem' higher resolution than the
    // normal
    // macro problem

    auto cpu_time = DSC_PROFILER.stopTiming("msfem.all");
    auto max_cpu_time = Dune::Fem::MPIManager::comm().max(cpu_time);
    DSC_LOG_INFO_0 << "Maximum total runtime of the program over all processes: " << max_cpu_time << "ms" << std::endl;
    DSC_PROFILER.outputTimings("profiler");

    // Compute the peak memory consumption of each processes
    int who = RUSAGE_SELF;
    struct rusage usage;
    getrusage(who, &usage);
    long peakMemConsumption = usage.ru_maxrss;
    // compute the maximum and mean peak memory consumption over all processes
    long maxPeakMemConsumption = Dune::Fem::MPIManager::comm().max(peakMemConsumption);
    long totalPeakMemConsumption = Dune::Fem::MPIManager::comm().sum(peakMemConsumption);
    long meanPeakMemConsumption = totalPeakMemConsumption / Dune::Fem::MPIManager::size();
    // write output on rank zero
    if (Dune::Fem::MPIManager::rank() == 0) {
      std::unique_ptr<boost::filesystem::ofstream> memoryConsFile(
          DSC::make_ofstream(DSC_CONFIG_GET("global.datadir", "data") + "/memory.csv"));
      *memoryConsFile << "global.maxPeakMemoryConsumption,global.meanPeakMemoryConsumption\n" << maxPeakMemConsumption
                      << "," << meanPeakMemConsumption << std::endl;
    }

    return 0;
  }
  catch (Dune::Exception& e) {
    std::cerr << e.what() << std::endl;
  }
  catch (const std::exception& ex) {
    std::cerr << "Caught std::exception: " << ex.what() << "\n";
  }
  catch (const std::string& ex) {
    std::cerr << "Caught string-type exception: " << ex << "\n";
  }
  catch (...) {
    std::cerr << "Exception of non-known type thrown!\n";
  }
  return 1;
} // main
