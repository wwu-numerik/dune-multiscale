message(AUTHOR_WARNING "TODO: Implement module test.")

include(DuneUtils)
include(DuneMPI)

find_package(METIS)
find_package(ParMETIS )
include(AddMETISFlags)
include(AddParMETISFlags)
include(CheckEmplace)

find_package(Eigen3 3.2.0)
find_package(TBB REQUIRED)

find_package(UMFPack REQUIRED)
include(AddUMFPackFlags)

find_package(FFTW REQUIRED)
include_directories(${FFTW_INCLUDES})

set(GRIDDIM 3 CACHE STRING "")
set(USE_ISTL_BACKEND 1 CACHE BOOLEAN "use dune-istl as the la-backend. Disable to use Eigen3 instead.")
set(USE_FEM_BACKEND 0 CACHE BOOLEAN "use dune-fem as the discretization-backend. Disable to use dune-pdelab instead.")

If    ("${CMAKE_BUILD_TYPE}" MATCHES "^REL")
  ADD_IF_SUPPORTED(CMAKE_CXX_FLAGS_RELEASE "-funroll-loops" "-m64" "-mfpmath=sse" "-falign-loops" "-mtune=native" "-march=native" 
  "-pipe" "-fomit-frame-pointer" "-O4" "-fno-alias" )
EndIf ("${CMAKE_BUILD_TYPE}" MATCHES "^REL")

add_definitions("-DSPGRID" "-DENABLE_PARMETIS=1" "-DGRIDDIM=${GRIDDIM}" -DMETISNAMEL )