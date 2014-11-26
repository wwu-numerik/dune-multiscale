// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include <config.h>
#include <thread>
#include <dune/multiscale/common/main_init.hh>
#include <dune/multiscale/msfem/algorithm.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/stuff/common/parallel/helper.hh>
#include <dune/stuff/common/profiler.hh>

#include <tbb/task_scheduler_init.h>

int main(int argc, char** argv) {
  using namespace Dune::Multiscale;
  try {  
    init(argc, argv);

    //!TODO include base in config
    DSC_PROFILER.startTiming("msfem.all");

    const std::string datadir = DSC_CONFIG_GET("global.datadir", "data/");

    // generate directories for data output
    DSC::testCreateDirectory(datadir);
    DSC_LOG_INFO_0 << boost::format("Data will be saved under: %s\nLogs will be saved under: %s/%s/ms.log.log\n") %
                          datadir % datadir % DSC_CONFIG_GET("logging.dir", "log");

    msfem_algorithm();

    auto comm = Dune::MPIHelper::getCollectiveCommunication();
    auto cpu_time = DSC_PROFILER.stopTiming("msfem.all", DSC_CONFIG_GET("global.output_walltime", false));
    auto max_cpu_time = comm.max(cpu_time);
    DSC_LOG_INFO_0 << "Maximum total runtime of the program over all processes: " << max_cpu_time << "ms" << std::endl;
    DSC_PROFILER.outputTimings("profiler");
    mem_usage();
  }
  catch (Dune::Exception& e) {
    return handle_exception(e);
  }
  catch (std::exception& s) {
    return handle_exception(s);
  }
  return 0;
} // main
