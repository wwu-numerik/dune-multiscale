if( NOT CMAKE_BUILD_TYPE )
	message( STATUS "No build type selected, setting default build type 'RELWITHDEBINFO'" )
	set( CMAKE_BUILD_TYPE "RELWITHDEBINFO" CACHE STRING "DEBUG/RELEASE")
endif( NOT CMAKE_BUILD_TYPE )

PROJECT(dune-multiscale CXX)

#set minimum cmake version
cmake_minimum_required(VERSION 2.8.6)

# mac-specific settings
IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  # Xcode-specific settings
  IF(CMAKE_GENERATOR MATCHES "Xcode")
    set(CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LANGUAGE_STANDARD "c++11")
    set(CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LIBRARY "libc++")  
  ENDIF(CMAKE_GENERATOR MATCHES "Xcode")
  LINK_DIRECTORIES( /opt/local/lib/ )
ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin") 

if(NOT (dune-common_DIR OR dune-common_ROOT OR
        "${CMAKE_PREFIX_PATH}" MATCHES ".*dune-common.*"))
    string(REPLACE  ${CMAKE_PROJECT_NAME} dune-common dune-common_DIR
      ${PROJECT_BINARY_DIR})
endif()

#find dune-common and set the module path
find_package(dune-common REQUIRED)
list(APPEND CMAKE_MODULE_PATH ${dune-common_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake/)

#include the dune macros
include(DuneMacros)

# start a dune project with information from dune.module
dune_project()
include(DuneUtils)
include(DuneMPI)

find_package(METIS)
find_package(ParMETIS REQUIRED)
include(cmake/modules/AddMETISFlags.cmake)
include(AddParMETISFlags)
include(FindTBB)
find_package(TBB )

foreach(_dep ${ALL_DEPENDENCIES})
  FILE( GLOB_RECURSE dune "${CMAKE_CURRENT_SOURCE_DIR}/../${_dep}/dune/*.hh" )
  set( DUNE_HEADERS ${DUNE_HEADERS} ${dune} )
endforeach(_dep ${ALL_DEPENDENCIES})

FILE( GLOB_RECURSE header "${CMAKE_CURRENT_SOURCE_DIR}/src/*.hh" )
FILE( GLOB_RECURSE dune "${CMAKE_CURRENT_SOURCE_DIR}/dune/*.hh" )

set( COMMON_HEADER ${header} ${dune} ${DUNE_HEADERS} )
set_source_files_properties( ${COMMON_HEADER} PROPERTIES HEADER_FILE_ONLY 1 )

#these should be split into all+hmm+msfem parts
set( COMMON_SOURCES 
        dune/multiscale/tools/misc/outputparameter.cc
        dune/multiscale/msfem/fem_solver.cc
        dune/multiscale/common/grid_creation.cc
        dune/multiscale/common/main_init.cc
        dune/multiscale/common/traits.cc
        dune/multiscale/fem/print_info.cc
        # error_calc
        dune/multiscale/common/error_calc.cc
        dune/multiscale/msfem/localsolution_proxy.cc
        dune/multiscale/msfem/localproblems/localgridsearch.cc
        # 
        dune/multiscale/msfem/proxygridview.cc
        dune/multiscale/msfem/localproblems/localgridlist.cc
    )

set( PROBLEM_SOURCES
    dune/multiscale/problems/elliptic/nine.cc
    dune/multiscale/problems/elliptic/spe10.cc
    dune/multiscale/problems/elliptic/tarbert.cc
    dune/multiscale/problems/selector.cc
    dune/multiscale/problems/base.cc
)

set( MSFEM_SOURCES
        dune/multiscale/msfem/algorithm.cc
        dune/multiscale/msfem/msfem_solver.cc
        dune/multiscale/msfem/coarse_scale_operator.cc
        
        dune/multiscale/msfem/localproblems/localoperator.cc
        dune/multiscale/msfem/localproblems/localproblemsolver.cc
        dune/multiscale/msfem/localproblems/localsolutionmanager.cc
        
        dune/multiscale/msfem/coarse_scale_assembler.cc
        dune/multiscale/msfem/coarse_rhs_functional.cc
    )

set( FEM_SOURCES
        dune/multiscale/fem/algorithm.cc
        src/elliptic_fem.cc
    )

If    ("${CMAKE_BUILD_TYPE}" MATCHES "^REL")
  ADD_IF_SUPPORTED(CMAKE_CXX_FLAGS_RELEASE "-funroll-loops" "-m64" "-mfpmath=sse" "-falign-loops" "-mtune=native" "-march=native" 
  "-pipe" "-fomit-frame-pointer" "-O3" "-fno-alias" )
EndIf ("${CMAKE_BUILD_TYPE}" MATCHES "^REL")

set(ENABLE_PETSC 0 CACHE BOOLEAN "")
set(GRIDDIM 2 CACHE STRING "")
set(USE_ISTL_BACKEND 1 CACHE BOOLEAN "use dune-istl as the la-backend. Disable to use Eigen3 instead.")
set(USE_FEM_BACKEND 0 CACHE BOOLEAN "use dune-fem as the discretization-backend. Disable to use dune-pdelab instead.")

add_definitions("-DSPGRID" "-DENABLE_PARMETIS=1" "-DGRIDDIM=${GRIDDIM}" ${PETSC_DEFINITIONS} -DMETISNAMEL )
ADD_LIBRARY(multiscale_common STATIC ${COMMON_SOURCES})
ADD_LIBRARY(multiscale_problem STATIC ${PROBLEM_SOURCES})
ADD_LIBRARY(multiscale_msfem STATIC ${MSFEM_SOURCES})

add_dune_mpi_flags(OBJECT multiscale_common)
add_dune_mpi_flags(OBJECT multiscale_problem)
add_dune_mpi_flags(OBJECT multiscale_msfem)

add_dune_metis_flags(multiscale_common multiscale_problem multiscale_msfem)
add_dune_parmetis_flags(multiscale_common multiscale_problem multiscale_msfem)

set( COMMON_LIBS multiscale_common multiscale_problem multiscale_common ${PETSC_LIBRARIES} ${DUNE_DEFAULT_LIBS})

include_directories( ${PETSC_INCLUDES} ${TBB_INCLUDE_DIRS} )
set(ALL_MULTISCALE_LIBS multiscale_msfem  ${COMMON_LIBS} )

ADD_EXECUTABLE(elliptic_msfem src/elliptic_msfem.cc ${COMMON_HEADER} )
TARGET_LINK_LIBRARIES(elliptic_msfem multiscale_msfem ${COMMON_LIBS} )
add_dune_mpi_flags(elliptic_msfem)

ADD_EXECUTABLE(elliptic_fem ${FEM_SOURCES} ${COMMON_HEADER} )
TARGET_LINK_LIBRARIES(elliptic_fem ${COMMON_LIBS} )
add_dune_mpi_flags(elliptic_fem)

ADD_EXECUTABLE(fail EXCLUDE_FROM_ALL src/fail.cc ${COMMON_HEADER})
# TARGET_LINK_LIBRARIES(fail ${COMMON_LIBS} )
add_dune_mpi_flags(fail)

add_dune_parmetis_flags(elliptic_msfem elliptic_fem fail)

ADD_EXECUTABLE(thread_vis EXCLUDE_FROM_ALL src/thread_vis.cc ${COMMON_HEADER} )
TARGET_LINK_LIBRARIES(thread_vis multiscale_msfem ${COMMON_LIBS} )
add_dune_mpi_flags(thread_vis)

add_dune_alugrid_flags(OBJECT multiscale_common)
add_dune_alugrid_flags(OBJECT multiscale_problem)
add_dune_alugrid_flags(OBJECT multiscale_msfem)
add_dune_alugrid_flags(elliptic_msfem)
add_dune_alugrid_flags(elliptic_fem)
add_dune_alugrid_flags(thread_vis)

# HEADERCHECK( ${header} ${dune} )
DEPENDENCYCHECK( ${header} ${dune} )

add_subdirectory(doc)
add_subdirectory(dune)

#disable most warnings from dependent modules
foreach(_mod ${ALL_DEPENDENCIES})
      dune_module_to_uppercase(_upper_case "${_mod}")
      if(${_mod}_INCLUDE_DIRS)
        foreach( _idir ${${_mod}_INCLUDE_DIRS} )
            add_definitions("-isystem ${_idir}")
        endforeach( _idir )
      endif(${_mod}_INCLUDE_DIRS)
endforeach(_mod DEPENDENCIES)
  
ADD_CPPCHECK(${CMAKE_CURRENT_SOURCE_DIR}/src/elliptic_msfem.cc 
             ${dune} ${header})
add_analyze(${COMMON_SOURCES} ${PROBLEM_SOURCES} ${MSFEM_SOURCES} 
            ${FEM_SOURCES} src/elliptic_msfem.cc )
add_format(${COMMON_SOURCES} ${PROBLEM_SOURCES} ${MSFEM_SOURCES} 
           ${FEM_SOURCES} src/elliptic_msfem.cc
           ${header} ${dune})

finalize_dune_project(GENERATE_CONFIG_H_CMAKE)

