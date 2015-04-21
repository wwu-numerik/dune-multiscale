#include <config.h>
#include "main_init.hh"

// for rusage
#include <sys/resource.h>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/exceptions.hh>

#if HAVE_DUNE_SPGRID
#include <dune/grid/spgrid/capabilities.hh>
#endif

#include <dune/stuff/common/configuration.hh>
#include <dune/stuff/common/logging.hh>
#include <dune/stuff/common/profiler.hh>
#include <dune/stuff/common/misc.hh>
#include <dune/stuff/common/parallel/threadmanager.hh>

#include <dune/multiscale/common/traits.hh>

void Dune::Multiscale::init(int argc, char** argv) {
#if DUNE_MULTISCALE_WITH_DUNE_FEM
  Dune::Fem::MPIManager::initialize(argc, argv);
#endif
  auto& helper = Dune::MPIHelper::instance(argc, argv);
  if (helper.size() > 1 && !(Dune::Capabilities::isParallel<Dune::Multiscale::CommonTraits::GridType>::v)) {
    DUNE_THROW(Dune::InvalidStateException, "mpi enabled + serial grid = bad idea");
  }
  DSC::Config().read_command_line(argc, argv);
  DSC::testCreateDirectory(DSC_CONFIG_GET("global.datadir", "data/"));

  // LOG_NONE = 1, LOG_ERROR = 2, LOG_INFO = 4,LOG_DEBUG = 8,LOG_CONSOLE = 16,LOG_FILE = 32
  // --> LOG_ERROR | LOG_INFO | LOG_DEBUG | LOG_CONSOLE | LOG_FILE = 62
  const bool useLogger = false;
  DSC::Logger().create(DSC_CONFIG_GETB("logging.level", 62, useLogger),
                       DSC_CONFIG_GETB("logging.file", std::string(argv[0]) + ".log", useLogger),
                       DSC_CONFIG_GETB("global.datadir", "data", useLogger),
                       DSC_CONFIG_GETB("logging.dir", "log" /*path below datadir*/, useLogger));
  DSC_CONFIG.set_record_defaults(true);
  DSC_PROFILER.setOutputdir(DSC_CONFIG_GET("global.datadir", "data"));
  DS::threadManager().set_max_threads(DSC_CONFIG_GET("threading.max_count", 4));
}

int Dune::Multiscale::handle_exception(const Dune::Exception& exp) {
  std::cerr << "Failed with Dune::Exception: " << exp.what();
  DSC_PROFILER.outputTimings("profiler");
  mem_usage();
  return Dune::Stuff::abort_all_mpi_processes();
}

int Dune::Multiscale::handle_exception(const std::exception& exp) {
  std::cerr << "Failed with std::exception: " << exp.what();
  DSC_PROFILER.outputTimings("profiler");
  mem_usage();
  return Dune::Stuff::abort_all_mpi_processes();
}

int Dune::Multiscale::handle_exception(const tbb::tbb_exception& exp) {
  std::cerr << "Failed with " << exp.name() << ": " << exp.what();
  DSC_PROFILER.outputTimings("profiler");
  mem_usage();
  return Dune::Stuff::abort_all_mpi_processes();
}

void Dune::Multiscale::mem_usage() {
  auto comm = Dune::MPIHelper::getCollectiveCommunication();
  // Compute the peak memory consumption of each processes
  int who = RUSAGE_SELF;
  struct rusage usage;
  getrusage(who, &usage);
  long peakMemConsumption = usage.ru_maxrss;
  // compute the maximum and mean peak memory consumption over all processes
  long maxPeakMemConsumption = comm.max(peakMemConsumption);
  long totalPeakMemConsumption = comm.sum(peakMemConsumption);
  long meanPeakMemConsumption = totalPeakMemConsumption / comm.size();
  // write output on rank zero
  if (comm.rank() == 0) {
    std::unique_ptr<boost::filesystem::ofstream> memoryConsFile(
        DSC::make_ofstream(std::string(DSC_CONFIG_GET("global.datadir", "data/")) + std::string("/memory.csv")));
    *memoryConsFile << "global.maxPeakMemoryConsumption,global.meanPeakMemoryConsumption\n" << maxPeakMemConsumption
                    << "," << meanPeakMemConsumption << std::endl;
  }
}


void Dune::Multiscale::dump_environment()
{
  std::unique_ptr<boost::filesystem::ofstream> of(DSC::make_ofstream(
                         std::string(DSC_CONFIG_GET("global.datadir", "data/")) + std::string("/env.txt")));
  DSC::dump_environment(*of);
}
