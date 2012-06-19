PROJECT(dune_multiscale)

#set minimum cmake version
cmake_minimum_required(VERSION 2.6)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake ${CMAKE_CURRENT_SOURCE_DIR}/../dune-stuff/cmake )
include(DuneUtils)
include(UseDoxygen)

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
	"YASPGRID" CACHE STRING
	"GRIDTYPE" )

SET( GRIDDIM
	"2" CACHE STRING
	"GRIDDIM" )

SET( ALUGRID_BASE_PATH
	"/opt/dune/modules/alugrid/alugrid-1.50" CACHE STRING
	"ALUGRID_BASE_PATH" )

SET( UGGRID_BASE_PATH
	"/share/dune/Modules/modules_x86_64/ug" CACHE STRING
	"UGGRID_BASE_PATH" )

SET_CONFIGHEADER_VARS()
CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/cmake/cmake_config.h.in ${CMAKE_CURRENT_BINARY_DIR}/cmake_config.h )

INCLUDE_DIR( ${dune_multiscale_BINARY_DIR} )
INCLUDE_DIR( ${dune_multiscale_SOURCE_DIR} )
INCLUDE_DIR( ${dune_multiscale_SOURCE_DIR}/.. )

IF( ${GRIDTYPE} STREQUAL "UGGRID" )
	ADD_CXX_FLAGS( -DENABLE_UG )
	INCLUDE_SYS_DIR( ${UGGRID_BASE_PATH}/include )
	LINK_DIRECTORIES( ${UGGRID_BASE_PATH}/lib/ )
	SET( UGLIBS "ugS3" "ugS2" "devS")
ENDIF( ${GRIDTYPE} STREQUAL "UGGRID" )

IF( ${ALUGRID_FOUND} )
	INCLUDE_SYS_DIR( ${ALUGRID_BASE_PATH}/include )
	INCLUDE_SYS_DIR( ${ALUGRID_BASE_PATH}/include/serial  )
	INCLUDE_SYS_DIR( ${ALUGRID_BASE_PATH}/include/duneinterface  )
	SET( ALUGRID_LIBS "alugrid" )
ENDIF( ${ALUGRID_FOUND} )	

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
	SET( ALULIBS "alugrid" )
	SET( ENABLE_ALUGRID ON )
	IF ( ENABLE_PARALELL )
		INCLUDE_SYS_DIR( ${ALUGRID_BASE_PATH}/include/parallel )
	ENDIF ( ENABLE_PARALELL )
	ADD_CXX_FLAGS( "-DENABLE_ALUGRID" )
ENDIF( ${GRIDTYPE} MATCHES "ALUGRID" )


CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/cmake/cmake_config.h.in ${CMAKE_CURRENT_BINARY_DIR}/cmake_config.h )

ADD_DUNE_MODULES( common istl grid fem geometry localfunctions subgrid )
ADD_MY_MODULES( stuff )
ADD_CXX_FLAGS( ${CUSTOM_FLAGS} )

LINK_DIRECTORIES( "/usr/lib" "${ALUGRID_BASE_PATH}/lib" )

FILE( GLOB_RECURSE header "src/*.hh" )
FILE( GLOB_RECURSE dune "dune/*.hh" )

set( COMMON_LIBS "dunefem" "dunegrid" "dunecommon" "dunegeometry" ${BLAS_LIB} ${ALUGRID_LIBS} ${UGLIBS} ${ParaLIBS} ${CCGNU_LIBRARIES} "boost_date_time" "boost_filesystem" "boost_system" )
set( COMMON_HEADER ${header} ${stokes} ${stuff} ${dune} ${DUNE_HEADERS} )
set_source_files_properties( ${COMMON_HEADER} PROPERTIES HEADER_FILE_ONLY 1 )

ADD_DEFINITIONS( ${MY_CXX_FLAGS} )

ADD_EXECUTABLE(elliptic_hmm src/elliptic_hmm.cc ${COMMON_HEADER} )
TARGET_LINK_LIBRARIES(elliptic_hmm ${COMMON_LIBS} )

ADD_EXECUTABLE(elliptic_msfem src/elliptic_msfem.cc ${COMMON_HEADER} )
TARGET_LINK_LIBRARIES(elliptic_msfem ${COMMON_LIBS} )

# HEADERCHECK( ${header} ${oseen} ${dune} )
# ADD_CPPCHECK( ${CMAKE_CURRENT_SOURCE_DIR}/src/elliptic_hmm.cc ${COMMON_HEADER} )

