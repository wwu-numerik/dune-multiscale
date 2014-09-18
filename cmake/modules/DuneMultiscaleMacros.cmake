message(AUTHOR_WARNING "TODO: Implement module test.")

find_package(PETSc COMPONENTS CXX )

find_package(METIS)
find_package(ParMETIS REQUIRED)
include(cmake/modules/AddMETISFlags.cmake)
include(AddParMETISFlags)
include(CheckEmplace)

find_package(Eigen3 3.2.0 REQUIRED)
find_package(TBB REQUIRED)