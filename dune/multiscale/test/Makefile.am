GRIDTYPE=YASPGRID
GRIDDIM=2

# tests where program to build and program to run are equal
NORMALTESTS = 

# list of tests to run
TESTS = $(NORMALTESTS)

# programs just to build when "make check" is used
check_PROGRAMS = $(NORMALTESTS)


LDFLAGS = -lboost_filesystem -lboost_system 
CXXFLAGS = $(BOOST_CPPFLAGS) $(DUNE_CPPFLAGS) \
	-DGRIDDIM=$(GRIDDIM) -D$(GRIDTYPE)


## distribution tarball
SOURCES = 

# gridcheck not used explicitly, we should still ship it :)
EXTRA_DIST = $(SOURCES)

CLEANFILES = *.gcda *.gcno semantic.cache simplex-testgrid*.dgf.* cube-testgrid*.dgf.* dgfparser.log

include $(top_srcdir)/am/global-rules
