// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_MULTISCALE_SRC_COMMON_HH
#define DUNE_MULTISCALE_SRC_COMMON_HH

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/exceptions.hh>

#if HAVE_DUNE_SPGRID
#include <dune/grid/spgrid/capabilities.hh>
#endif

#include <dune/stuff/common/configuration.hh>
#include <dune/stuff/common/logging.hh>
#include <dune/stuff/common/profiler.hh>
#include <dune/stuff/common/parallel/threadmanager.hh>

#include <dune/multiscale/common/traits.hh>

namespace Dune {
namespace Multiscale {

//! setup code common to fem/msfem/hmm
void init(int argc, char** argv) {
  namespace DSC = Dune::Stuff::Common;
#if DUNE_MULTISCALE_WITH_DUNE_FEM
  Dune::Fem::MPIManager::initialize(argc, argv);
#endif
  auto& helper = Dune::MPIHelper::instance(argc, argv);
  if (helper.size() > 1 && !(Dune::Capabilities::isParallel<Dune::Multiscale::CommonTraits::GridType>::v)) {
    DUNE_THROW(Dune::InvalidStateException, "mpi enabled + serial grid = bad idea");
  }
  DSC::Config().read_command_line(argc, argv);

  // LOG_NONE = 1, LOG_ERROR = 2, LOG_INFO = 4,LOG_DEBUG = 8,LOG_CONSOLE = 16,LOG_FILE = 32
  // --> LOG_ERROR | LOG_INFO | LOG_DEBUG | LOG_CONSOLE | LOG_FILE = 62
  const bool useLogger = false;
  DSC::Logger().create(DSC_CONFIG_GETB("logging.level", 62, useLogger),
                       DSC_CONFIG_GETB("logging.file", std::string(argv[0]) + ".log", useLogger),
                       DSC_CONFIG_GETB("global.datadir", "data", useLogger),
                       DSC_CONFIG_GETB("logging.dir", "log" /*path below datadir*/, useLogger));
  DSC_CONFIG.set_record_defaults(true);
  DSC_PROFILER.setOutputdir(DSC_CONFIG_GET("global.datadir", "data"));
  DS::ThreadManager::set_max_threads(DSC_CONFIG_GET("threading.max_count", 1));
} // init

} // namespace Dune {
} // namespace Multiscale {

#endif // DUNE_MULTISCALE_SRC_COMMON_HH
