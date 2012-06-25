#ifndef DUNE_MULTISCALE_SRC_COMMON_HH
#define DUNE_MULTISCALE_SRC_COMMON_HH

#ifdef HAVE_CMAKE_CONFIG
 #include "cmake_config.h"
#elif defined (HAVE_CONFIG_H)
 #include <config.h>
#endif // ifdef HAVE_CMAKE_CONFIG

// polynomial order of discrete space
#define POLORDER 1

#ifndef USE_GRAPE
 #define USE_GRAPE HAVE_GRAPE
#endif

#define USE_TWISTFREE_MAPPER
#define VERBOSE false

#include <iostream>
#include <sstream>

#include <stdio.h>
#include <stdlib.h>
// -----------------------------

#include <dune/common/mpihelper.hh> // An initializer of MPI
#include <dune/common/exceptions.hh> // We use exceptions
#include <dune/common/unused.hh>

#if HAVE_GRAPE
 #include <dune/grid/io/visual/grapedatadisplay.hh>
#endif

// to display data with ParaView:
#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/fem/io/file/dataoutput.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/io/file/datawriter.hh>
#include <dune/fem/gridpart/common/gridpart.hh>
#include <dune/multiscale/grids/periodicgridpart/periodicgridpart.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/space/lagrangespace.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/operator/matrix/spmatrix.hh>
#include <dune/fem/space/common/adaptmanager.hh>
#include <dune/fem/solver/inverseoperators.hh>

#include <dune/stuff/configcontainer.hh>
#include <dune/stuff/parametercontainer.hh>
#include <dune/stuff/debug.hh>
#include <dune/stuff/misc.hh>
#include <dune/stuff/logging.hh>

#if ENABLE_MPI
        typedef Dune::CollectiveCommunication< MPI_Comm > CollectiveCommunication;
#else
        typedef Dune::CollectiveCommunication< double > CollectiveCommunication;
#endif

CollectiveCommunication init( int argc, char** argv )
{
    Dune::MPIManager::initialize(argc, argv);
    Stuff::Config().readCommandLine( argc, argv );

    // LOG_NONE = 1, LOG_ERR = 2, LOG_INFO = 4,LOG_DEBUG = 8,LOG_CONSOLE = 16,LOG_FILE = 32
    //--> LOG_ERR | LOG_INFO | LOG_DEBUG | LOG_CONSOLE | LOG_FILE = 62
    const bool useLogger = false;
    Logger().Create( Stuff::Config().get( "logging.level",  62,                             useLogger ),
                     Stuff::Config().get( "logging.file",   std::string(argv[0]) + ".log",  useLogger ),
                     Stuff::Config().get( "global.datadir",   "data",                         useLogger ),
                     Stuff::Config().get( "logging.dir",    ""/*path below datadir*/,       useLogger )
                    );

    return CollectiveCommunication();//( Dune::MPIManager::helper().getCommunicator() );
}

#endif // DUNE_MULTISCALE_SRC_COMMON_HH
