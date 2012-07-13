
SUBDIRS =

GRIDTYPE=ALBERTAGRID
GRIDDIM=2
PROBLEM_NAMES=Easy Toy One Two Three Four Five Six Seven Eight Nine Ten
HMM_PROBLEM=One
MSFEM_PROBLEM=Ten

noinst_PROGRAMS = elliptic_hmm elliptic_msfem 
SOURCE_PREFIX=elliptic_hmm elliptic_msfem
elliptic_hmm_SOURCES = elliptic_hmm.cc
elliptic_msfem_SOURCES = elliptic_msfem.cc

common_CPPFLAGS = $(ALL_PKG_CPPFLAGS) $(AM_CPPFLAGS) \
	$(DUNEMPICPPFLAGS) \
	$(UG_CPPFLAGS) \
	$(AMIRAMESH_CPPFLAGS) \
	$(ALBERTA_CPPFLAGS) \
	$(ALUGRID_CPPFLAGS)  \
	-DGRIDDIM=$(GRIDDIM) -DGRIDTYPE=$(GRIDTYPE) -DWORLDDIM=$(GRIDDIM)
	
elliptic_hmm_CPPFLAGS=$(common_CPPFLAGS) -DPROBLEM_NAME=$(HMM_PROBLEM)
elliptic_msfem_CPPFLAGS = $(common_CPPFLAGS) -DPROBLEM_NAME=$(MSFEM_PROBLEM)

# The libraries have to be given in reverse order (most basic libraries
# last).  Also, due to some misunderstanding, a lot of libraries include the
# -L option in LDFLAGS instead of LIBS -- so we have to include the LDFLAGS
# here as well.
common_LDADD = \
	$(DUNE_LDFLAGS) $(DUNE_LIBS) \
	$(ALUGRID_LDFLAGS) $(ALUGRID_LIBS) \
	$(ALBERTA_LDFLAGS) $(ALBERTA_LIBS) \
	$(AMIRAMESH_LDFLAGS) $(AMIRAMESH_LIBS) \
	$(UG_LDFLAGS) $(UG_LIBS) -lboost_filesystem -lboost_system \
	$(DUNEMPILIBS)	\
	$(LDADD)
common_LDFLAGS = $(AM_LDFLAGS) \
	$(DUNEMPILDFLAGS) \
	$(UG_LDFLAGS) \
	$(AMIRAMESH_LDFLAGS) \
	$(ALBERTA_LDFLAGS) \
	$(ALUGRID_LDFLAGS) \
	$(DUNE_LDFLAGS)

elliptic_hmm_LDADD=$(common_LDADD)
elliptic_msfem_LDADD=$(common_LDADD)
elliptic_hmm_LDFLAGS=$(common_LDFLAGS)
elliptic_msfem_LDFLAGS=$(common_LDFLAGS)

.PHONY: all_problems
all_problems: 
	@for problem in $(PROBLEM_NAMES) ; do \
		$(MAKE) clean; \
		$(MAKE) PROBLEM_NAME=$$problem elliptic_hmm; \
		$(MAKE) PROBLEM_NAME=$$problem elliptic_msfem; \
	done

TESTS=all_problems

# don't follow the full GNU-standard
# we need automake 1.5
AUTOMAKE_OPTIONS = foreign 1.5

# pass most important options when "make distcheck" is used
DISTCHECK_CONFIGURE_FLAGS = --with-dune-grid=$(DUNE_GRID_ROOT) --with-dune-fem=$(DUNE_FEM_ROOT) --with-dune-common=$(DUNE_COMMON_ROOT) --with-ag-common=$(AG_COMMON_ROOT) --with-dune-istl=$(DUNE_ISTL_ROOT) --with-dune-localfunctions=$(DUNE_LOCALFUNCTIONS_ROOT)  CXX="$(CXX)" CC="$(CC)"

include $(top_srcdir)/am/global-rules

