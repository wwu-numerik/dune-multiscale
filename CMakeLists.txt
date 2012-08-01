PROJECT(dune_multiscale)

#set minimum cmake version
cmake_minimum_required(VERSION 2.6)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake ${CMAKE_CURRENT_SOURCE_DIR}/../dune-stuff/cmake )
include(DuneUtils)

#----------------------------------------------------------------------------------------------------
# General Settings
#----------------------------------------------------------------------------------------------------
SET( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true )
SET( CMAKE_FIND_LIBRARY_SUFFIXES ".so" ".lib" ".la" ".a")
SET( CMAKE_COLOR_MAKEFILE ON)

EXECUTE_PROCESS(
	COMMAND ${CMAKE_SOURCE_DIR}/cmake/getRevision.sh
	OUTPUT_VARIABLE COMMIT
	ERROR_VARIABLE shell_error
	OUTPUT_STRIP_TRAILING_WHITESPACE
)

#----------------------------------------------------------------------------------------------------
# Options, that can be changed be the User in order to customise dune_multiscale
#----------------------------------------------------------------------------------------------------
SET( BLAS_LIB
	"blas" CACHE STRING
	"" )

SET( USE_ISTL
	0 CACHE STRING
	"Use fem-wrapped istl matrix" )
		
SET( GRIDTYPE
	"ALBERTAGRID" CACHE STRING
	"GRIDTYPE" )

SET( GRIDDIM
	"2" CACHE STRING
	"GRIDDIM" )

SET( ALUGRID_BASE_PATH
	"/opt/dune/modules/alugrid/alugrid-1.50" CACHE STRING
	"ALUGRID_BASE_PATH" )
SET( ALBERTA_GRID_BASE_PATH
	"/opt/dune/modules/alberta/alberta-3.0" CACHE STRING
	"ALUGRID_BASE_PATH" )
SET( UGGRID_BASE_PATH
	"/share/dune/Modules/modules_x86_64/ug" CACHE STRING
	"UGGRID_BASE_PATH" )

SET(PROBLEM_NAMES "Easy" "Toy" "One" "Two" "Three" "Four" "Five" "Six" "Seven" "Eight" "Nine" "Ten" INTERNAL )

SET( HMM_PROBLEM
	"Easy" CACHE STRING
	"HMM_PROBLEM NAME" )

SET( MSFEM_PROBLEM
	"Ten" CACHE STRING
	"MSFEM_PROBLEM NAME" )



SET_CONFIGHEADER_VARS()
CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/cmake/cmake_config.h.in ${CMAKE_CURRENT_BINARY_DIR}/cmake_config.h )

INCLUDE_DIR( ${dune_multiscale_BINARY_DIR} )
INCLUDE_DIR( ${dune_multiscale_SOURCE_DIR} )
INCLUDE_DIR( ${dune_multiscale_SOURCE_DIR}/.. )

IF( ${GRIDTYPE} STREQUAL "UGGRID" )
	ADD_CXX_FLAGS( -DENABLE_UG )
	INCLUDE_SYS_DIR( ${UGGRID_BASE_PATH}/include )
	LINK_DIRECTORIES( ${UGGRID_BASE_PATH}/lib/ )
	SET( GRIDLIBS "ugS3" "ugS2" "devS")
ENDIF( ${GRIDTYPE} STREQUAL "UGGRID" )


SET( ENABLE_ALUGRID OFF )
IF( ${GRIDTYPE} MATCHES "ALUGRID" )
#	IF( ${ALUGRID_FOUND} )
#		MESSAGE( FATAL ) 
#	ENDIF( ${ALUGRID_FOUND} )	
	INCLUDE_SYS_DIR( 	${GRAPE_PATH} 
						${ALUGRID_BASE_PATH}/include 
						${ALUGRID_BASE_PATH}/include/serial
						${ALUGRID_BASE_PATH}/include/duneinterface )
	LINK_DIRECTORIES( "${ALUGRID_BASE_PATH}/lib"  )
	SET( GRIDLIBS "alugrid" )
	SET( ENABLE_ALUGRID ON )
	IF ( ENABLE_PARALELL )
		INCLUDE_SYS_DIR( ${ALUGRID_BASE_PATH}/include/parallel )
	ENDIF ( ENABLE_PARALELL )
	ADD_CXX_FLAGS( "-DENABLE_ALUGRID" )
ENDIF( ${GRIDTYPE} MATCHES "ALUGRID" )

SET( ENABLE_ALBERTA OFF )
IF( ${GRIDTYPE} MATCHES "ALBERTAGRID" )
	INCLUDE_SYS_DIR( 	${ALBERTA_GRID_BASE_PATH}/include 
						${ALBERTA_GRID_BASE_PATH}/include/alberta)
	LINK_DIRECTORIES( "${ALBERTA_GRID_BASE_PATH}/lib"  )
# 	SET( GRIDLIBS alberta_utilities_debug alberta_2d_debug dunealbertagrid_2d dunegrid)
	SET( GRIDLIBS alberta_utilities alberta_2d dunealbertagrid_2d dunegrid)
	SET( ENABLE_ALBERTA ON )
	ADD_CXX_FLAGS( -DENABLE_ALBERTA=1 -DDEBUG )
ENDIF( ${GRIDTYPE} MATCHES "ALBERTAGRID" )

CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/cmake/cmake_config.h.in ${CMAKE_CURRENT_BINARY_DIR}/cmake_config.h )

ADD_DUNE_MODULES( common istl grid fem geometry localfunctions subgrid stuff )
ADD_CXX_FLAGS( ${CUSTOM_FLAGS} -DGRIDTYPE=${GRIDTYPE} -DHAVE_CMAKE_CONFIG -DHAVE_CONFIG_H)

LINK_DIRECTORIES( "/usr/lib" "${ALUGRID_BASE_PATH}/lib" )

FILE( GLOB_RECURSE header "src/*.hh" )
FILE( GLOB_RECURSE dune "dune/*.hh" )

set( COMMON_LIBS "dunefem" "dunegrid" "dunecommon" "dunegeometry" ${BLAS_LIB} ${GRIDLIBS} ${ParaLIBS} ${CCGNU_LIBRARIES} "boost_date_time" "boost_filesystem" "boost_system" )
set( COMMON_HEADER ${header} ${stokes} ${stuff} ${dune} ${DUNE_HEADERS} )
set_source_files_properties( ${COMMON_HEADER} PROPERTIES HEADER_FILE_ONLY 1 )

ADD_DEFINITIONS( ${MY_CXX_FLAGS} )

ADD_EXECUTABLE(elliptic_hmm src/elliptic_hmm.cc ${COMMON_HEADER} )
TARGET_LINK_LIBRARIES(elliptic_hmm ${COMMON_LIBS} )
set_target_properties(elliptic_hmm PROPERTIES COMPILE_FLAGS -DPROBLEM_NAME=${HMM_PROBLEM})


ADD_EXECUTABLE(elliptic_msfem src/elliptic_msfem.cc ${COMMON_HEADER} )
TARGET_LINK_LIBRARIES(elliptic_msfem ${COMMON_LIBS} )
set_target_properties(elliptic_msfem PROPERTIES COMPILE_FLAGS -DPROBLEM_NAME=${MSFEM_PROBLEM})

ADD_EXECUTABLE(msfem-1d src/msfem-1d.cc ${COMMON_HEADER} )
TARGET_LINK_LIBRARIES(msfem-1d ${COMMON_LIBS} )

# HEADERCHECK( ${header} ${oseen} ${dune} )
ADD_CPPCHECK( ${CMAKE_CURRENT_SOURCE_DIR}/src/elliptic_hmm.cc ${CMAKE_CURRENT_SOURCE_DIR}/src/elliptic_msfem.cc ${COMMON_HEADER} )

add_subdirectory(doc)