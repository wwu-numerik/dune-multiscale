message(AUTHOR_WARNING "TODO: Implement module test.")

find_package(PETSc COMPONENTS CXX )

find_package(METIS)
find_package(ParMETIS)
include(cmake/modules/AddMETISFlags.cmake)
include(AddParMETISFlags)