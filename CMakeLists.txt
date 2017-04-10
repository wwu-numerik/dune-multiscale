PROJECT(dune-multiscale CXX)

# local environment
set( ENV{PATH} "${CMAKE_CURRENT_SOURCE_DIR}/../local/bin:$ENV{PATH}" )
set( ENV{LD_LIBRARY_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/../local/lib:$ENV{LD_LIBRARY_PATH}" )
set( ENV{PKG_CONFIG_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/../local/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}" )

#set minimum cmake version
cmake_minimum_required(VERSION 3.1)
set( CMAKE_COLOR_MAKEFILE ON)

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
list(APPEND CMAKE_MODULE_PATH ${dune-common_MODULE_PATH}
"${PROJECT_SOURCE_DIR}/cmake/modules")

#include the dune macros
include(DuneMacros)

# start a dune project with information from dune.module
dune_project()
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules )

foreach(_dep ${ALL_DEPENDENCIES})
  FILE( GLOB_RECURSE dune "${CMAKE_CURRENT_SOURCE_DIR}/../${_dep}/dune/*.hh" )
  set( DUNE_HEADERS ${DUNE_HEADERS} ${dune} )
endforeach(_dep ${ALL_DEPENDENCIES})

FILE( GLOB_RECURSE header "${CMAKE_CURRENT_SOURCE_DIR}/src/*.hh" )
FILE( GLOB_RECURSE dune "${CMAKE_CURRENT_SOURCE_DIR}/dune/*.hh" )

set( COMMON_HEADER ${header} ${dune} ${DUNE_HEADERS} )
set_source_files_properties( ${COMMON_HEADER} PROPERTIES HEADER_FILE_ONLY 1 )

dune_enable_all_packages(
   INCLUDE_DIRS ${dune-multiscale}/dune
   COMPILE_DEFINITIONS
   COMPILE_OPTIONS
   MODULE_LIBRARIES multiscale_common multiscale_problem multiscale_cgfem multiscale_msfem
)

# macros from dune-xt
include(DuneUtils)

#these should be split into all+hmm+msfem parts
set( COMMON_SOURCES
        dune/multiscale/tools/misc/outputparameter.cc
        dune/multiscale/msfem/fem_solver.cc
        dune/multiscale/common/grid_creation.cc
        dune/multiscale/common/main_init.cc
        dune/multiscale/common/traits.cc
        dune/multiscale/common/df_io.cc
        dune/multiscale/fem/print_info.cc
        # error_calc
        dune/multiscale/common/error_calc.cc
        dune/multiscale/common/heterogenous.cc
        dune/multiscale/msfem/localsolution_proxy.cc
        dune/multiscale/msfem/localproblems/localgridsearch.cc
        #
        dune/multiscale/msfem/proxygridview.cc
        dune/multiscale/msfem/localproblems/localgridlist.cc
    )

set( PROBLEM_SOURCES
    dune/multiscale/problems/synthetic.cc
    dune/multiscale/problems/er2007.cc
    dune/multiscale/problems/spe10.cc
    dune/multiscale/problems/tarbert.cc
    dune/multiscale/problems/selector.cc
    dune/multiscale/problems/base.cc
    dune/multiscale/problems/random.cc
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

set( CGFEM_SOURCES
        dune/multiscale/fem/algorithm.cc
    )

dune_library_add_sources(multiscale_common SOURCES ${COMMON_SOURCES})
dune_library_add_sources(multiscale_problem SOURCES ${PROBLEM_SOURCES})
dune_library_add_sources(multiscale_msfem SOURCES ${MSFEM_SOURCES})
dune_library_add_sources(multiscale_cgfem SOURCES ${CGFEM_SOURCES})

ADD_EXECUTABLE(elliptic_msfem src/msfem.cc ${COMMON_HEADER} )
ADD_EXECUTABLE(elliptic_fem src/cgfem.cc ${CGFEM_SOURCES} ${COMMON_HEADER} )
ADD_EXECUTABLE(fail EXCLUDE_FROM_ALL src/fail.cc ${COMMON_HEADER})
ADD_EXECUTABLE(minimal EXCLUDE_FROM_ALL src/minimal.cc ${COMMON_HEADER})
ADD_EXECUTABLE(thread_vis EXCLUDE_FROM_ALL src/thread_vis.cc ${COMMON_HEADER} )

DEPENDENCYCHECK( ${header} ${dune} )

add_subdirectory(doc)
add_subdirectory(dune)

#disable most warnings from dependent modules, except gdt,stuff
list(APPEND module_exceptions "dune-stuff" "dune-xt-common" "dune-xt-la" "dune-xt-grid" "dune-xt-functions" "dune-gdt")
foreach(_mod ${ALL_DEPENDENCIES})
  list (FIND module_exceptions "${_mod}" _index)
  if (_index LESS 0)
    dune_module_to_uppercase(_upper_case "${_mod}")
    if(${_mod}_INCLUDE_DIRS)
      foreach( _idir ${${_mod}_INCLUDE_DIRS} )
        add_definitions("-isystem ${_idir}")
      endforeach( _idir )
    endif(${_mod}_INCLUDE_DIRS)
  endif()
endforeach(_mod DEPENDENCIES)

add_format(${CMAKE_CURRENT_SOURCE_DIR})
add_tidy(${CMAKE_CURRENT_SOURCE_DIR})

if (NOT ${SPE10MODEL1DATA} STREQUAL "SPE10MODEL1DATA-NOTFOUND")
  set (SPE10MODEL1DATA_TARGET_FILENAME "${CMAKE_CURRENT_BINARY_DIR}/spe10_permeability.dat")
  if (NOT EXISTS ${SPE10MODEL1DATA_TARGET_FILENAME})
    execute_process(COMMAND ln -s "${SPE10MODEL1DATA}" "${SPE10MODEL1DATA_TARGET_FILENAME}")
  endif(NOT EXISTS ${SPE10MODEL1DATA_TARGET_FILENAME})
else(NOT ${SPE10MODEL1DATA} STREQUAL "SPE10MODEL1DATA-NOTFOUND")
  message(STATUS "Downloading Spe10 model1 data file:" ${SPE10MODEL1DATA})
  execute_process(COMMAND wget -nc http://wwwmath.uni-muenster.de/num/ohlberger/research/software/perm_case1.zip
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  execute_process(COMMAND unzip -o perm_case1.zip WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  set (SPE10MODEL1DATA_TARGET_FILENAME "${CMAKE_CURRENT_BINARY_DIR}/spe10_permeability.dat")
  execute_process(COMMAND mv perm_case1.dat "${SPE10MODEL1DATA_TARGET_FILENAME}")
endif (NOT ${SPE10MODEL1DATA} STREQUAL "SPE10MODEL1DATA-NOTFOUND")

finalize_dune_project(GENERATE_CONFIG_H_CMAKE)

