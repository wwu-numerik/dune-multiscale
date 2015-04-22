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

# guess dune-common build dir
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
        dune/multiscale/common/heterogenous.cc
        dune/multiscale/msfem/localsolution_proxy.cc
        dune/multiscale/msfem/localproblems/localgridsearch.cc
        # 
        dune/multiscale/msfem/proxygridview.cc
        dune/multiscale/msfem/localproblems/localgridlist.cc
        dune/multiscale/tools/discretefunctionwriter.cc
    )

set( PROBLEM_SOURCES
    dune/multiscale/problems/synthetic.cc
    dune/multiscale/problems/spe10.cc
    dune/multiscale/problems/tarbert.cc
    dune/multiscale/problems/random.cc
    dune/multiscale/problems/selector.cc
    dune/multiscale/problems/base.cc
)
if(HAVE_RANDOM_PROBLEM)
    set( PROBLEM_SOURCES ${PROBLEM_SOURCES} )
endif(HAVE_RANDOM_PROBLEM)

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

dune_add_library(multiscale_common ${COMMON_SOURCES} ADD_LIBS ${DUNE_LIBS} ${COMMON_LIBS})
dune_add_library(multiscale_problem ${PROBLEM_SOURCES} ADD_LIBS ${DUNE_LIBS} ${COMMON_LIBS})
dune_add_library(multiscale_msfem ${MSFEM_SOURCES} ADD_LIBS ${DUNE_LIBS} ${COMMON_LIBS})
message(STATUS MULTI ${MULTISCALE_LIBS})
ADD_EXECUTABLE(elliptic_msfem src/elliptic_msfem.cc ${COMMON_HEADER} )
TARGET_LINK_LIBRARIES(elliptic_msfem multiscale_msfem ${MULTISCALE_LIBS} )

ADD_EXECUTABLE(elliptic_fem ${FEM_SOURCES} ${COMMON_HEADER} )
TARGET_LINK_LIBRARIES(elliptic_fem ${MULTISCALE_LIBS} )

ADD_EXECUTABLE(fail EXCLUDE_FROM_ALL src/fail.cc ${COMMON_HEADER})
# TARGET_LINK_LIBRARIES(fail ${MULTISCALE_LIBS} )

ADD_EXECUTABLE(thread_vis EXCLUDE_FROM_ALL src/thread_vis.cc ${COMMON_HEADER} )
TARGET_LINK_LIBRARIES(thread_vis multiscale_msfem ${MULTISCALE_LIBS} )

foreach(target multiscale_common multiscale_problem multiscale_msfem elliptic_msfem elliptic_fem thread_vis)
  add_dune_metis_flags(${target})
  add_dune_parmetis_flags(${target})
  add_dune_superlu_flags(${target})
  add_dune_umfpack_flags(${target})
  add_dune_tbb_flags(${target})
  add_dune_alugrid_flags(${target})
  add_dune_mpi_flags(${target})
endforeach()

# HEADERCHECK( ${header} ${dune} )
DEPENDENCYCHECK( ${header} ${dune} )

add_subdirectory(doc)
add_subdirectory(dune)

#disable most warnings from dependent modules, except gdt,stuff
list(APPEND module_exceptions "dune-stuff" "dune-gdt")
foreach(_mod ${ALL_DEPENDENCIES})
  list (FIND module_exceptions "${_mod}" _index)
  if (_index LESS 0)
    dune_module_to_uppercase(_upper_case "${_mod}")
    if(${_mod}_INCLUDE_DIRS)
      foreach( _idir ${${_mod}_INCLUDE_DIRS} )
        add_definitions("-isystem ${_idir}")
      endforeach( _idir )
    endif(${_mod}_INCLUDE_DIRS)
  endif(_index GREATER -1)
endforeach(_mod DEPENDENCIES)
  
ADD_CPPCHECK(${CMAKE_CURRENT_SOURCE_DIR}/src/elliptic_msfem.cc 
            src/thread_vis.cc ${dune} ${header})
add_analyze(${COMMON_SOURCES} ${PROBLEM_SOURCES} ${MSFEM_SOURCES} 
            ${FEM_SOURCES} src/elliptic_msfem.cc src/thread_vis.cc )
add_format(${COMMON_SOURCES} ${PROBLEM_SOURCES} ${MSFEM_SOURCES} 
           ${FEM_SOURCES} src/elliptic_msfem.cc src/thread_vis.cc
           ${header} ${dune})

finalize_dune_project(GENERATE_CONFIG_H_CMAKE)

