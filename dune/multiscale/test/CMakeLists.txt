BEGIN_TESTCASES()


foreach( test ${testnames} )
    add_dune_mpi_flags(${test})
    add_dune_alugrid_flags(${test})
    TARGET_LINK_LIBRARIES(${test} multiscale_msfem ${COMMON_LIBS} )
endforeach( test ${testnames} )

END_TESTCASES()
