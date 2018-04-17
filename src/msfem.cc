// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include <config.h>
#include <thread>
#include <dune/multiscale/common/main_init.hh>
#include <dune/multiscale/msfem/algorithm.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/xt/common/parallel/helper.hh>
#include <dune/xt/common/timings.hh>
#include <dune/xt/common/filesystem.hh>
#include <dune/xt/common/configuration.hh>

#include <tbb/task_scheduler_init.h>

int main(int argc, char** argv)
{
  using namespace Dune::Multiscale;
  // try {
    init(argc, argv);
    const size_t max_threads = DXTC_CONFIG_GET("threading.max_count", 1);
    tbb::task_scheduler_init sched_init(max_threads);
    Dune::XT::Common::threadManager().set_max_threads(max_threads);

    //! TODO include base in config
    DXTC_TIMINGS.start("msfem.all");

    const std::string datadir = DXTC_CONFIG_GET("global.datadir", "data/");

    // generate directories for data output
    Dune::XT::Common::test_create_directory(datadir);
    MS_LOG_INFO_0 << boost::format("Data will be saved under: %s\nLogs will be saved under: %s/%s/ms.log.log\n")
                         % datadir % datadir % DXTC_CONFIG_GET("logging.dir", "log");

    msfem_algorithm();

    auto comm = Dune::MPIHelper::getCollectiveCommunication();
    auto cpu_time = DXTC_TIMINGS.stop("msfem.all");
    auto max_cpu_time = comm.max(cpu_time);
    MS_LOG_INFO_0 << "Maximum total runtime of the program over all processes: " << max_cpu_time << "ms" << std::endl;
    DXTC_TIMINGS.output_per_rank("profiler");
    mem_usage();
    dump_environment();
  // } catch (Dune::Exception& e) {
  //   return handle_exception(e);
  // } catch (std::exception& s) {
  //   return handle_exception(s);
  // }
  return 0;
} // main
