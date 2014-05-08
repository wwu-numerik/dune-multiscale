#include <config.h>

// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include <dune/multiscale/common/main_init.hh>
#include <dune/multiscale/common/error_container.hh>
#include <dune/multiscale/msfem/algorithm.hh>
#include <dune/multiscale/problems/selector.hh>
#include <dune/multiscale/tools/discretefunctionwriter.hh>

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

    algorithm();

    auto cpu_time = DSC_PROFILER.stopTiming("msfem.all", DSC_CONFIG_GET("global.output_walltime", false));
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
