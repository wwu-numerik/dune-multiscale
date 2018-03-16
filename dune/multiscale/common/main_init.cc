#include <config.h>
#include "main_init.hh"

// for rusage
#include <sys/resource.h>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/exceptions.hh>

#if HAVE_DUNE_SPGRID
#include <dune/grid/spgrid/capabilities.hh>
#endif

#include <dune/xt/common/configuration.hh>
#include <dune/xt/common/logging.hh>
#include <dune/xt/common/signals.hh>
#include <dune/xt/common/timings.hh>
#include <dune/xt/common/misc.hh>
#include <dune/xt/common/parallel/threadmanager.hh>
#include <dune/xt/common/filesystem.hh>

#include <dune/multiscale/common/traits.hh>

#include <tbb/tbb_exception.h>

static void handle_sigterm(int)
{
  DXTC_TIMINGS.stop();
  DXTC_TIMINGS.output_per_rank("profiler");
  std::exit(5);
}

void Dune::Multiscale::init(int argc, char** argv)
{
#if DUNE_MULTISCALE_WITH_DUNE_FEM
  Dune::Fem::MPIManager::initialize(argc, argv);
#endif
  auto&& helper = Dune::MPIHelper::instance(argc, argv);
  //  if (helper.size() > 1 && !(Dune::Capabilities::isParallel<Dune::Multiscale::CommonTraits::GridType>::v)) {
  //    DUNE_THROW(Dune::InvalidStateException, "mpi enabled + serial grid = bad idea");
  //  }
  // makes inserting '\n' not flush
  std::cout.sync_with_stdio(false);

  Dune::XT::Common::Config().read_command_line(argc, argv);
  Dune::XT::Common::Config().read_command_line(argc, argv);
  Dune::XT::Common::test_create_directory(DXTC_CONFIG_GET("global.datadir", "data/"));

  // LOG_NONE = 1, LOG_ERROR = 2, LOG_INFO = 4,LOG_DEBUG = 8,LOG_CONSOLE = 16,LOG_FILE = 32
  // --> LOG_ERROR | LOG_INFO | LOG_DEBUG | LOG_CONSOLE | LOG_FILE = 62
  Dune::XT::Common::Logger().create(DXTC_CONFIG_GET("logging.level", 62),
                                    DXTC_CONFIG_GET("logging.file", std::string(argv[0]) + ".log"),
                                    DXTC_CONFIG_GET("global.datadir", "data"),
                                    DXTC_CONFIG_GET("logging.dir", "log" /*path below datadir*/));
  DXTC_TIMINGS.set_outputdir(DXTC_CONFIG_GET("global.datadir", "data"));

  Dune::XT::Common::install_signal_handler(SIGTERM, handle_sigterm);
#ifdef MS_TIMED_LOGGER
  Dune::XT::Common::TimedLogger().create(20, // max info level
                                         20, // max debug level
                                         true, // warnings are enabled (the default)
                                         true, // colors are enabled (the default)
                                         "white", // info color (the default)
                                         "lightgrey", // debug color (the default)
                                         "red" // warning color (the default)
                                         );
#endif
}

int Dune::Multiscale::handle_exception(const Dune::Exception& exp)
{
  std::cerr << "Failed with Dune::Exception: " << exp.what();
  DXTC_TIMINGS.output_per_rank("profiler");
  mem_usage();
  return Dune::XT::abort_all_mpi_processes();
}

int Dune::Multiscale::handle_exception(const std::exception& exp)
{
  std::cerr << "Failed with std::exception: " << exp.what();
  DXTC_TIMINGS.output_per_rank("profiler");
  mem_usage();
  return Dune::XT::abort_all_mpi_processes();
}

int Dune::Multiscale::handle_exception(const tbb::tbb_exception& exp)
{
  std::cerr << "Failed with tbb::exception" << exp.name() << ": " << exp.what();
  DXTC_TIMINGS.output_per_rank("profiler");
  mem_usage();
  return Dune::XT::abort_all_mpi_processes();
}

void Dune::Multiscale::mem_usage()
{
  auto comm = Dune::MPIHelper::getCollectiveCommunication();
  // Compute the peak memory consumption of each processes
  int who = RUSAGE_SELF;
  struct rusage usage;
  getrusage(who, &usage);
  long peakMemConsumption = usage.ru_maxrss;
  // compute the maximum and mean peak memory consumption over all processes
  long maxPeakMemConsumption = comm.max(peakMemConsumption);
  const long totalPeakMemConsumption = comm.sum(peakMemConsumption);
  const long meanPeakMemConsumption = totalPeakMemConsumption / comm.size();
  // write output on rank zero
  if (comm.rank() == 0) {
    std::unique_ptr<boost::filesystem::ofstream> memoryConsFile(Dune::XT::Common::make_ofstream(
        std::string(DXTC_CONFIG_GET("global.datadir", "data/")) + std::string("/memory.csv")));
    *memoryConsFile << "global.maxPeakMemoryConsumption,global.meanPeakMemoryConsumption\n"
                    << maxPeakMemConsumption << "," << meanPeakMemConsumption << std::endl;
  }
}

void Dune::Multiscale::dump_environment()
{
  std::unique_ptr<boost::filesystem::ofstream> of(Dune::XT::Common::make_ofstream(
      std::string(DXTC_CONFIG_GET("global.datadir", "data/")) + std::string("/env.txt")));
  Dune::XT::Common::dump_environment(*of);
}
