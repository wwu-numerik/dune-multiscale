
set(TESTS  compare_hmm compare_msfem compare_lodm)
# set(LDFLAGS  -lboost_timer -lboost_chrono -lboost_filesystem -lboost_system ${BOOST_LDFLAGS})
# set(CXXFLAGS  ${BOOST_CPPFLAGS} ${DUNE_CPPFLAGS}
#     -DGRIDDIM=${GRIDDIM} -D${GRIDTYPE} ${ALL_PKG_CPPFLAGS})
add_executable(compare_hmm compare_hmm.cc)
add_executable(compare_msfem compare_msfem.cc)
add_executable(compare_lodm compare_lodm.cc)

foreach(i ${TESTS})
  add_test(${i} ${i})
endforeach(i ${TESTS})
# We do not want want to build the tests during make all
# but just build them on demand
add_directory_test_target(_test_target)
add_dependencies(${_test_target} ${TESTS} ${COMPILE_XFAIL_TESTS})