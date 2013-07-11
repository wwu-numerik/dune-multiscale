ENABLE_TESTING()

BEGIN_TESTCASES()

foreach( algo hmm msfem lodm )
  target_link_libraries( test_compare_${algo} ${ALL_MULTISCALE_LIBS}  pthread )
endforeach( algo )


END_TESTCASES()