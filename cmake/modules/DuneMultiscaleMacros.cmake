message(AUTHOR_WARNING "TODO: Implement module test.")

include(DuneUtils)
include(DuneMPI)

find_package(METIS)
find_package(ParMETIS )
include(cmake/modules/AddMETISFlags.cmake)
include(AddParMETISFlags)
include(CheckEmplace)

find_package(Eigen3 3.2.0 REQUIRED)
find_package(TBB REQUIRED)

find_package(UMFPack REQUIRED)
include(AddUMFPackFlags)
