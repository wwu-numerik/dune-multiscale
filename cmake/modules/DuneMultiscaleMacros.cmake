set(DUNE_GRID_GRIDTYPE_SELECTOR ON)
include(DuneUtils)
include(DuneMPI)

find_package(METIS)
find_package(ParMETIS )
include(AddMETISFlags)
include(AddParMETISFlags)
include(CheckEmplace)

find_package(UMFPack REQUIRED)
find_package(SuiteSparse REQUIRED)
include_directories( ${SUITESPARSE_INCLUDE_DIRS} )
set(DUNE_UMFPACK_LIBRARIES ${UMFPACK_LIBRARY} ${CHOLMOD_LIBRARY} ${COLAMD_LIBRARY} ${AMD_LIBRARY} ${SUITESPARSE_CONFIG_LIBRARY} )

find_package(FFTW REQUIRED)
set(HAVE_RANDOM_PROBLEM ${HAVE_FFTW})

set( MULTISCALE_LIBS multiscale_common multiscale_cgfem multiscale_msfem multiscale_problem multiscale_common ${DUNE_DEFAULT_LIBS} ${COMMON_LIBS} ${DUNE_UMFPACK_LIBRARIES} )

dune_define_gridtype(GRID_CONFIG_H_BOTTOM
                     GRIDTYPE YASPGRID_OFFSET
                     ASSERTION "GRIDDIM == WORLDDIM"
                     DUNETYPE "Dune::YaspGrid< dimgrid, Dune::EquidistantOffsetCoordinates<double, dimgrid> >"
                     HEADERS "dune/grid/yaspgrid.hh" "dune/grid/io/file/dgfparser/dgfyasp.hh")

set(MSFEM_MACRO_GRIDDIM 3 CACHE STRING "Make sure to use the same settings in dune-mlmc")
set(MSFEM_MACRO_GRIDTYPE YASPGRID_OFFSET CACHE STRING "AnyOf SPGRID_ISOTROPIC SPGRID_ANISOTROPIC SPGRID_BISECTION YASPGRID_OFFSET")
set(USE_ISTL_BACKEND 1 CACHE BOOLEAN "use dune-istl as the la-backend. Disable to use Eigen3 instead.")
set(USE_FEM_BACKEND 0 CACHE BOOLEAN "use dune-fem as the discretization-backend. Disable to use dune-pdelab instead.")

dune_register_package_flags(COMPILE_DEFINITIONS "ENABLE_PARMETIS=1;GRIDDIM=${MSFEM_MACRO_GRIDDIM};${MSFEM_MACRO_GRIDTYPE};USE_ISTL_BACKEND=${USE_ISTL_BACKEND};USE_FEM_BACKEND=${USE_FEM_BACKEND}")
