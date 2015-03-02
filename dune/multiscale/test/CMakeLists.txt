BEGIN_TESTCASES()


foreach( test ${testnames} )
  add_dune_mpi_flags(${test})
  add_dune_alugrid_flags(${test})
  add_dune_metis_flags(${test})
  add_dune_parmetis_flags(${test})
  add_dune_superlu_flags(${test})
  add_dune_umfpack_flags(${test})
  add_dune_alugrid_flags(${test})
  add_dune_mpi_flags(${test})
  target_link_libraries(${test} multiscale_msfem ${MULTISCALE_LIBS} ${COMMON_LIBS} )
endforeach( test ${testnames} )

END_TESTCASES()
