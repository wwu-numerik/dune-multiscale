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

# dune_define_gridtype(GRID_CONFIG_H_BOTTOM
#                      GRIDTYPE YASPGRID_OFFSET
#                      ASSERTION "GRIDDIM == WORLDDIM"
#                      DUNETYPE "Dune::YaspGrid< dimgrid, Dune::EquidistantOffsetCoordinates<double, world_dim> >"
#                      HEADERS "dune/grid/yaspgrid.hh" "dune/grid/io/file/dgfparser/dgfyasp.hh")
                     
