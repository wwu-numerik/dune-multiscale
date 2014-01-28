
# add_executable(compare_hmm compare_hmm.cc)

BEGIN_TESTCASES()

TARGET_LINK_LIBRARIES(test_compare_msfem multiscale_msfem ${COMMON_LIBS} )
TARGET_LINK_LIBRARIES(test_compare_hmm multiscale_hmm ${COMMON_LIBS} )

foreach( test ${testnames} )
    add_dune_mpi_flags(${test})
    add_dune_alugrid_flags(${test})
endforeach( test ${testnames} )

END_TESTCASES()