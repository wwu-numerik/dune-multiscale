message(AUTHOR_WARNING "TODO: Implement module test.")

if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  find_package(PETSc COMPONENTS CXX REQUIRED)
else(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  find_package(PETSc REQUIRED)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
