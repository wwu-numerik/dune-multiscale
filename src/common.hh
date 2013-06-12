// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_MULTISCALE_SRC_COMMON_HH
#define DUNE_MULTISCALE_SRC_COMMON_HH

#ifdef HAVE_CMAKE_CONFIG
 #include "cmake_config.h"
#elif defined (HAVE_CONFIG_H)
 #include <config.h>
#endif // ifdef HAVE_CMAKE_CONFIG

// polynomial order of discrete space
#ifndef POLORDER
  #define POLORDER 1
#endif

#ifndef USE_GRAPE
 #define USE_GRAPE HAVE_GRAPE
#endif

#define USE_TWISTFREE_MAPPER

#include <iostream>
#include <sstream>

#include <stdio.h>
#include <stdlib.h>
// -----------------------------

#include <dune/common/mpihelper.hh> // An initializer of MPI
#include <dune/common/exceptions.hh> // We use exceptions
#include <dune/common/unused.hh>

// to display data with ParaView:
#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/fem/io/file/dataoutput.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/io/file/datawriter.hh>
#include <dune/fem/gridpart/common/gridpart.hh>
#include <dune/fem/gridpart/periodicgridpart/periodicgridpart.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/space/lagrangespace.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/operator/matrix/spmatrix.hh>
#include <dune/fem/space/common/adaptmanager.hh>
#include <dune/fem/solver/inverseoperators.hh>
#include <dune/fem/misc/mpimanager.hh>

#include <dune/stuff/common/parameter/configcontainer.hh>
#include <dune/stuff/common/debug.hh>
#include <dune/stuff/common/misc.hh>
#include <dune/stuff/common/logging.hh>
#include <dune/stuff/common/profiler.hh>
#include <dune/stuff/aliases.hh>

#include <dune/multiscale/common/traits.hh>

void init(int argc, char** argv) {
  namespace DSC = Dune::Stuff::Common;
  Dune::MPIManager::initialize(argc, argv);
  if (Dune::MPIManager::size() > 1
      && !(Dune::Capabilities::isParallel<Dune::Multiscale::CommonTraits::GridType>::v))
  {
    DUNE_THROW(Dune::InvalidStateException, "mpi enabled + serial grid = bad idea");
  }
  DSC::Config().readCommandLine(argc, argv);

  // LOG_NONE = 1, LOG_ERROR = 2, LOG_INFO = 4,LOG_DEBUG = 8,LOG_CONSOLE = 16,LOG_FILE = 32
  // --> LOG_ERROR | LOG_INFO | LOG_DEBUG | LOG_CONSOLE | LOG_FILE = 62
  const bool useLogger = false;
  DSC::Logger().create(DSC_CONFIG_GETB("logging.level", 62, useLogger),
                  DSC_CONFIG_GETB("logging.file", std::string(argv[0]) + ".log", useLogger),
                  DSC_CONFIG_GETB("global.datadir", "data", useLogger),
                  DSC_CONFIG_GETB("logging.dir", "log" /*path below datadir*/, useLogger)
                  );
  DSC_CONFIG.setRecordDefaults(true);
  DSC_PROFILER.setOutputdir(DSC_CONFIG_GET("global.datadir", "data"));
} // init

#endif // DUNE_MULTISCALE_SRC_COMMON_HH
