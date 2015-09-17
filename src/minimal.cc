// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include <config.h>
#include <dune/multiscale/common/main_init.hh>
#include <dune/multiscale/fem/algorithm.hh>
#include <dune/stuff/common/parallel/helper.hh>
#include <dune/stuff/common/profiler.hh>
#include <thread>

#include <tbb/task_scheduler_init.h>

int main(int argc, char** argv) {
  using namespace Dune::Multiscale;
  try {
#if DUNE_MULTISCALE_WITH_DUNE_FEM
//  Dune::Fem::MPIManager::initialize(argc, argv);
#endif
  auto& helper = Dune::MPIHelper::instance(argc, argv);
//  if (helper.size() > 1 && !(Dune::Capabilities::isParallel<Dune::Multiscale::CommonTraits::GridType>::v)) {
//    DUNE_THROW(Dune::InvalidStateException, "mpi enabled + serial grid = bad idea");
//  }
//  DSC::Config().read_command_line(argc, argv);
//  DSC::testCreateDirectory(DSC_CONFIG_GET("global.datadir", "data/"));

//  // LOG_NONE = 1, LOG_ERROR = 2, LOG_INFO = 4,LOG_DEBUG = 8,LOG_CONSOLE = 16,LOG_FILE = 32
//  // --> LOG_ERROR | LOG_INFO | LOG_DEBUG | LOG_CONSOLE | LOG_FILE = 62
//  DSC::Logger().create(DSC_CONFIG_GET("logging.level", 62),
//                       DSC_CONFIG_GET("logging.file", std::string(argv[0]) + ".log"),
//                       DSC_CONFIG_GET("global.datadir", "data"),
//                       DSC_CONFIG_GET("logging.dir", "log" /*path below datadir*/));
//  DSC_PROFILER.setOutputdir(DSC_CONFIG_GET("global.datadir", "data"));
//  DS::threadManager().set_max_threads(DSC_CONFIG_GET("threading.max_count",1));

//    DSC_PROFILER.startTiming("total_cpu");

//    cgfem_algorithm();

//    const auto cpu_time =
//        DSC_PROFILER.stopTiming("total_cpu") / 1000.f;
//    DSC_LOG_INFO_0 << "Total runtime of the program: " << cpu_time << "s" << std::endl;
//    DSC_PROFILER.outputTimings("profiler");
//    mem_usage();
//    dump_environment();
  } catch (Dune::Exception& e) {
    return handle_exception(e);
  } catch (std::exception& s) {
    return handle_exception(s);
  }
  return 0;
} // main
