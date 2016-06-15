message(AUTHOR_WARNING "TODO: Implement module test.")

include(DuneUtils)
include(DuneMPI)

find_package(METIS)
find_package(ParMETIS )
include(AddMETISFlags)
include(AddParMETISFlags)
include(CheckEmplace)

find_package(UMFPack)
find_package(SuiteSparse)
include_directories( ${SUITESPARSE_INCLUDE_DIRS} )
set(DUNE_UMFPACK_LIBRARIES ${UMFPACK_LIBRARY} ${CHOLMOD_LIBRARY} ${COLAMD_LIBRARY} ${AMD_LIBRARY} ${SUITESPARSE_CONFIG_LIBRARY} ) 

find_package(FFTW)
set(HAVE_RANDOM_PROBLEM ${HAVE_FFTW})

set( MULTISCALE_LIBS multiscale_common multiscale_cgfem multiscale_msfem multiscale_problem multiscale_common ${DUNE_DEFAULT_LIBS} ${COMMON_LIBS} ${DUNE_UMFPACK_LIBRARIES} )

set(GRIDDIM 3 CACHE STRING "")
set(USE_ISTL_BACKEND 1 CACHE BOOLEAN "use dune-istl as the la-backend. Disable to use Eigen3 instead.")
set(USE_FEM_BACKEND 0 CACHE BOOLEAN "use dune-fem as the discretization-backend. Disable to use dune-pdelab instead.")

add_definitions("-DENABLE_PARMETIS=1" "-DGRIDDIM=${GRIDDIM}" "-DMETISNAMEL" )
