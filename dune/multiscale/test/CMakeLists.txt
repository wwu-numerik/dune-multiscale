ENABLE_TESTING()

BEGIN_TESTCASES()

foreach( algo hmm msfem lodm )
    set( targetName test_${algo} )
    add_executable(${targetName} EXCLUDE_FROM_ALL ${dune_multiscale_SOURCE_DIR}/src/elliptic_${algo}.cc ${COMMON_HEADER} )
    target_link_libraries( ${targetName} ${ALL_MULTISCALE_LIBS} )
    add_test( NAME compile_${targetName} COMMAND ${CMAKE_MAKE_PROGRAM} ${targetName} )
    list(APPEND testnames compile_${targetName} )
endforeach( algo )


target_link_libraries( test_compare_msfem ${ALL_MULTISCALE_LIBS}  pthread )

END_TESTCASES()