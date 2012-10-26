
GRIDTYPE=ALUGRID_SIMPLEX
GRIDDIM=2

# tests where program to build and program to run are equal
NORMALTESTS = test_hmm_easy test_hmm_toy test_hmm_one test_hmm_two test_hmm_three test_hmm_four test_hmm_five test_hmm_six test_hmm_seven test_hmm_eight test_hmm_nine test_hmm_ten test_msfem_easy test_msfem_toy test_msfem_one test_msfem_two test_msfem_three test_msfem_four test_msfem_five test_msfem_six test_msfem_seven test_msfem_eight test_msfem_nine test_msfem_ten 

# list of tests to run
TESTS = $(NORMALTESTS)

# programs just to build when "make check" is used
check_PROGRAMS = $(NORMALTESTS)


LDFLAGS = -lboost_timer -lboost_chrono -lboost_filesystem -lboost_system $(BOOST_LDFLAGS)
CXXFLAGS = $(BOOST_CPPFLAGS) $(DUNE_CPPFLAGS) \
    -DGRIDDIM=$(GRIDDIM) -D$(GRIDTYPE) $(ALL_PKG_CPPFLAGS)


test_hmm_easy_SOURCES = $(top_srcdir)/src/elliptic_hmm.cc
test_hmm_easy_CPPFLAGS = -DPROBLEM_NAME=Easy

test_hmm_toy_SOURCES = $(top_srcdir)/src/elliptic_hmm.cc
test_hmm_toy_CPPFLAGS = -DPROBLEM_NAME=Toy

test_hmm_one_SOURCES = $(top_srcdir)/src/elliptic_hmm.cc
test_hmm_one_CPPFLAGS = -DPROBLEM_NAME=One

test_hmm_two_SOURCES = $(top_srcdir)/src/elliptic_hmm.cc
test_hmm_two_CPPFLAGS = -DPROBLEM_NAME=Two

test_hmm_three_SOURCES = $(top_srcdir)/src/elliptic_hmm.cc
test_hmm_three_CPPFLAGS = -DPROBLEM_NAME=Three

test_hmm_four_SOURCES = $(top_srcdir)/src/elliptic_hmm.cc
test_hmm_four_CPPFLAGS = -DPROBLEM_NAME=Four

test_hmm_five_SOURCES = $(top_srcdir)/src/elliptic_hmm.cc
test_hmm_five_CPPFLAGS = -DPROBLEM_NAME=Five

test_hmm_six_SOURCES = $(top_srcdir)/src/elliptic_hmm.cc
test_hmm_six_CPPFLAGS = -DPROBLEM_NAME=Six

test_hmm_seven_SOURCES = $(top_srcdir)/src/elliptic_hmm.cc
test_hmm_seven_CPPFLAGS = -DPROBLEM_NAME=Seven

test_hmm_eight_SOURCES = $(top_srcdir)/src/elliptic_hmm.cc
test_hmm_eight_CPPFLAGS = -DPROBLEM_NAME=Eight

test_hmm_nine_SOURCES = $(top_srcdir)/src/elliptic_hmm.cc
test_hmm_nine_CPPFLAGS = -DPROBLEM_NAME=Nine

test_hmm_ten_SOURCES = $(top_srcdir)/src/elliptic_hmm.cc
test_hmm_ten_CPPFLAGS = -DPROBLEM_NAME=Ten

test_msfem_easy_SOURCES = $(top_srcdir)/src/elliptic_msfem.cc
test_msfem_easy_CPPFLAGS = -DPROBLEM_NAME=Easy

test_msfem_toy_SOURCES = $(top_srcdir)/src/elliptic_msfem.cc
test_msfem_toy_CPPFLAGS = -DPROBLEM_NAME=Toy

test_msfem_one_SOURCES = $(top_srcdir)/src/elliptic_msfem.cc
test_msfem_one_CPPFLAGS = -DPROBLEM_NAME=One

test_msfem_two_SOURCES = $(top_srcdir)/src/elliptic_msfem.cc
test_msfem_two_CPPFLAGS = -DPROBLEM_NAME=Two

test_msfem_three_SOURCES = $(top_srcdir)/src/elliptic_msfem.cc
test_msfem_three_CPPFLAGS = -DPROBLEM_NAME=Three

test_msfem_four_SOURCES = $(top_srcdir)/src/elliptic_msfem.cc
test_msfem_four_CPPFLAGS = -DPROBLEM_NAME=Four

test_msfem_five_SOURCES = $(top_srcdir)/src/elliptic_msfem.cc
test_msfem_five_CPPFLAGS = -DPROBLEM_NAME=Five

test_msfem_six_SOURCES = $(top_srcdir)/src/elliptic_msfem.cc
test_msfem_six_CPPFLAGS = -DPROBLEM_NAME=Six

test_msfem_seven_SOURCES = $(top_srcdir)/src/elliptic_msfem.cc
test_msfem_seven_CPPFLAGS = -DPROBLEM_NAME=Seven

test_msfem_eight_SOURCES = $(top_srcdir)/src/elliptic_msfem.cc
test_msfem_eight_CPPFLAGS = -DPROBLEM_NAME=Eight

test_msfem_nine_SOURCES = $(top_srcdir)/src/elliptic_msfem.cc
test_msfem_nine_CPPFLAGS = -DPROBLEM_NAME=Nine

test_msfem_ten_SOURCES = $(top_srcdir)/src/elliptic_msfem.cc
test_msfem_ten_CPPFLAGS = -DPROBLEM_NAME=Ten
    

## distribution tarball
SOURCES = $(top_srcdir)/src/elliptic_hmm.cc $(top_srcdir)/src/elliptic_msfem.cc 

# gridcheck not used explicitly, we should still ship it :)
EXTRA_DIST = $(SOURCES)

CLEANFILES = *.gcda *.gcno semantic.cache simplex-testgrid*.dgf.* cube-testgrid*.dgf.* dgfparser.log

include $(top_srcdir)/am/global-rules
