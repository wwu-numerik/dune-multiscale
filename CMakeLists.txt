if( NOT CMAKE_BUILD_TYPE )
	message( STATUS "No build type selected, setting default build type 'RELWITHDEBINFO'" )
	set( CMAKE_BUILD_TYPE "RELWITHDEBINFO" CACHE STRING "DEBUG/RELEASE")
endif( NOT CMAKE_BUILD_TYPE )

PROJECT(dune-multiscale CXX)

#set minimum cmake version
cmake_minimum_required(VERSION 2.8.6)

#find dune-common and set the module path
find_package(dune-common)
list(APPEND CMAKE_MODULE_PATH ${dune-common_MODULE_PATH}
     "${CMAKE_SOURCE_DIR}/cmake/modules" ${CMAKE_CURRENT_SOURCE_DIR}/cmake ${CMAKE_CURRENT_SOURCE_DIR}/../dune-stuff/cmake )
     # make sure our own modules are found

#include the dune macros
include(DuneMacros)

# start a dune project with information from dune.module
dune_project()
include(DuneUtils)

CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/cmake/cmake_config.h.in ${CMAKE_CURRENT_BINARY_DIR}/cmake_config.h )

IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
#	INCLUDE_SYS_DIR( opt/local/include )
#	INCLUDE_SYS_DIR( opt/local/include/openmpi )
  LINK_DIRECTORIES( /opt/local/lib/ )
ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin") 

FILE( GLOB_RECURSE header "${CMAKE_CURRENT_SOURCE_DIR}/src/*.hh" )
FILE( GLOB_RECURSE dune "${CMAKE_CURRENT_SOURCE_DIR}/dune/*.hh" )

set( COMMON_HEADER ${header} ${stokes} ${stuff} ${dune} ${DUNE_HEADERS} )
set_source_files_properties( ${COMMON_HEADER} PROPERTIES HEADER_FILE_ONLY 1 )

#these should be split into all+hmm+msfem parts
set( COMMON_SOURCES 
        dune/multiscale/tools/misc/outputparameter.cc
        dune/multiscale/hmm/cell_problem_numbering.cc
        dune/multiscale/hmm/cell_problem_solver.cc
        dune/multiscale/hmm/discrete_cell_operator.cc
        dune/multiscale/msfem/fem_solver.cc
        dune/multiscale/msfem/elliptic_msfem_matrix_assembler.cc
        dune/multiscale/msfem/localproblems/subgrid-list.cc
        dune/multiscale/msfem/msfem_grid_specifier.cc
        dune/multiscale/msfem/localproblems/localoperator.cc
        dune/multiscale/msfem/localproblems/localproblemsolver.cc
        dune/multiscale/common/error_calc.cc
        dune/multiscale/msfem/localproblems/weighted-clement-operator.cc
        dune/multiscale/msfem/localproblems/localsolutionmanager.cc
    )

set( PROBLEM_SOURCES
	dune/multiscale/problems/elliptic/eight.cc
	dune/multiscale/problems/elliptic/example.cc
	dune/multiscale/problems/elliptic/eleven.cc
	dune/multiscale/problems/elliptic/five.cc
	dune/multiscale/problems/elliptic/four.cc
	dune/multiscale/problems/elliptic/nine.cc
	dune/multiscale/problems/elliptic/one.cc
	dune/multiscale/problems/elliptic/seven.cc
	dune/multiscale/problems/elliptic/six.cc
	dune/multiscale/problems/elliptic/ten.cc
	dune/multiscale/problems/elliptic/thirteen.cc
	dune/multiscale/problems/elliptic/three.cc
	dune/multiscale/problems/elliptic/toy.cc
	dune/multiscale/problems/elliptic/twelve.cc
	dune/multiscale/problems/elliptic/two.cc
	dune/multiscale/problems/elliptic/spe10.cc
	dune/multiscale/problems/selector.cc
	dune/multiscale/problems/base.cc
    )
    
set( HMM_SOURCES  
        dune/multiscale/hmm/algorithm.cc
        dune/multiscale/hmm/algorithm_error.cc
        dune/multiscale/hmm/algorithm_step.cc
        dune/multiscale/hmm/elliptic_hmm_matrix_assembler.cc
    )

set( MSFEM_SOURCES
        dune/multiscale/msfem/algorithm.cc
        dune/multiscale/msfem/msfem_solver.cc
    )

set( HOM_FEM_SOURCES
        src/elliptic_hom_fem.cc
        dune/multiscale/fem/algorithm.cc
    )
    
set( FEM_SOURCES
        src/elliptic_fem.cc
        dune/multiscale/fem/algorithm.cc
    )
    
set( LODM_SOURCES 
        dune/multiscale/msfem/rigorous.cc
        dune/multiscale/msfem/rigorous_msfem_solver.cc
    )
    
If    ("${CMAKE_BUILD_TYPE}" MATCHES "^REL")
  ADD_IF_SUPPORTED(CMAKE_CXX_FLAGS_RELEASE "-funroll-loops" "-m64" "-mfpmath=sse" "-falign-loops" "-mtune=native" "-march=native" 
  "-pipe" "-fomit-frame-pointer" "-flto" "-O4" "-fwhole-program" "-ftree-vectorize" "-fno-alias")
EndIf ("${CMAKE_BUILD_TYPE}" MATCHES "^REL")

add_definitions("-DSPGRID" "-DGRIDDIM=2")
ADD_LIBRARY(multiscale_common STATIC ${COMMON_SOURCES})
ADD_LIBRARY(multiscale_problem STATIC ${PROBLEM_SOURCES})
ADD_LIBRARY(multiscale_lodm STATIC ${LODM_SOURCES})
ADD_LIBRARY(multiscale_msfem STATIC ${MSFEM_SOURCES})
ADD_LIBRARY(multiscale_hmm STATIC ${HMM_SOURCES})

set( COMMON_LIBS multiscale_problem multiscale_common "dunestuff" "dunefem" ${DUNE_DEFAULT_LIBS} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

set(ALL_MULTISCALE_LIBS multiscale_lodm multiscale_msfem multiscale_hmm ${COMMON_LIBS} )

ADD_EXECUTABLE(elliptic_hmm src/elliptic_hmm.cc ${COMMON_HEADER} )
TARGET_LINK_LIBRARIES(elliptic_hmm multiscale_hmm ${COMMON_LIBS} )

ADD_EXECUTABLE(elliptic_msfem src/elliptic_msfem.cc ${MSFEM_SOURCES} ${COMMON_HEADER} )
TARGET_LINK_LIBRARIES(elliptic_msfem multiscale_msfem ${COMMON_LIBS} )

ADD_EXECUTABLE(elliptic_hom_fem ${HOM_FEM_SOURCES} ${COMMON_HEADER} )
TARGET_LINK_LIBRARIES(elliptic_hom_fem ${COMMON_LIBS} )

ADD_EXECUTABLE(elliptic_fem ${FEM_SOURCES} ${COMMON_HEADER} )
TARGET_LINK_LIBRARIES(elliptic_fem ${COMMON_LIBS} )

ADD_EXECUTABLE(elliptic_lodm src/elliptic_lodm.cc  ${COMMON_HEADER} )
TARGET_LINK_LIBRARIES(elliptic_lodm multiscale_lodm ${COMMON_LIBS} )

ADD_EXECUTABLE(msfem-1d src/msfem-1d.cc ${COMMON_HEADER} )
TARGET_LINK_LIBRARIES(msfem-1d ${COMMON_LIBS} )


HEADERCHECK( ${header} ${dune} )
DEPENDENCYCHECK( ${header} ${dune} )

add_subdirectory(doc)
add_subdirectory(dune)

ADD_CPPCHECK( ${CMAKE_CURRENT_SOURCE_DIR}/src/elliptic_hmm.cc ${CMAKE_CURRENT_SOURCE_DIR}/src/elliptic_msfem.cc ${dune} ${header})
finalize_dune_project(GENERATE_CONFIG_H_CMAKE)

