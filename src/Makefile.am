
SUBDIRS =

#GRIDTYPE=ALBERTAGRID
GRIDTYPE=SPGRID
GRIDDIM=2
PROBLEM_NAMES= Toy One Two Three Four Five Six Seven Eight Nine Ten
FEM_PROBLEM=Nine
HMM_PROBLEM=Three
MSFEM_PROBLEM=Nine
MSFEM_RIGOROUS_PROBLEM=Nine
DEFAULT_PROBLEM=One

noinst_PROGRAMS = elliptic_hmm elliptic_msfem elliptic_fem elliptic_hom_fem elliptic_rigorous_msfem
SOURCE_PREFIX=elliptic_hmm elliptic_msfem elliptic_fem elliptic_hom_fem elliptic_rigorous_msfem
elliptic_hmm_SOURCES = elliptic_hmm.cc
elliptic_msfem_SOURCES = elliptic_msfem.cc
elliptic_rigorous_msfem_SOURCES = elliptic_rigorous_msfem.cc
elliptic_fem_SOURCES = elliptic_fem.cc
elliptic_hom_fem_SOURCES = elliptic_hom_fem.cc

CPPFLAGS = $(ALL_PKG_CPPFLAGS) $(AM_CPPFLAGS) \
	$(DUNEMPICPPFLAGS) \
	$(UG_CPPFLAGS) \
	$(AMIRAMESH_CPPFLAGS) \
	$(ALBERTA_CPPFLAGS) $(BOOST_CPPFLAGS) $(DUNE_CPPFLAGS) \
	$(ALUGRID_CPPFLAGS)  \
	-DGRIDDIM=$(GRIDDIM) -DGRIDTYPE=$(GRIDTYPE) -DWORLDDIM=$(GRIDDIM)
	
elliptic_hmm_CPPFLAGS=$(common_CPPFLAGS) -DPROBLEM_NAME=$(HMM_PROBLEM)
elliptic_msfem_CPPFLAGS = $(common_CPPFLAGS) -DPROBLEM_NAME=$(MSFEM_PROBLEM)
elliptic_rigorous_msfem_CPPFLAGS = $(common_CPPFLAGS) -DPROBLEM_NAME=$(MSFEM_RIGOROUS_PROBLEM)
elliptic_fem_CPPFLAGS=$(common_CPPFLAGS) -DPROBLEM_NAME=$(FEM_PROBLEM)
elliptic_hom_fem_CPPFLAGS = $(common_CPPFLAGS) -DPROBLEM_NAME=$(DEFAULT_PROBLEM)

# The libraries have to be given in reverse order (most basic libraries
# last).  Also, due to some misunderstanding, a lot of libraries include the
# -L option in LDFLAGS instead of LIBS -- so we have to include the LDFLAGS
# here as well.
LDADD = \
	$(DUNE_LDFLAGS) $(DUNE_LIBS) \
	$(ALUGRID_LDFLAGS) $(ALUGRID_LIBS) \
	$(ALBERTA_LDFLAGS) $(ALBERTA_LIBS) \
	$(AMIRAMESH_LDFLAGS) $(AMIRAMESH_LIBS) \
	$(UG_LDFLAGS) $(UG_LIBS) $(BOOST_LDFLAGS) \
	$(BOOST_FILESYSTEM_LIB) $(BOOST_SYSTEM_LIB) $(BOOST_TIMER_LIB) $(BOOST_CHRONO_LIB) \
	-lboost_filesystem-libc++ -lboost_system-libc++ -lboost_timer-libc++ -lboost_chrono-libc++ \
	$(DUNEMPILIBS)	\
	$(DUNE_LDADD)
LDFLAGS = $(AM_LDFLAGS) \
	$(DUNEMPILDFLAGS) \
	$(UG_LDFLAGS) \
	$(AMIRAMESH_LDFLAGS) \
	$(BOOST_LDFLAGS) \
	$(ALBERTA_LDFLAGS) \
	$(ALUGRID_LDFLAGS) \
	$(DUNE_LDFLAGS)

# don't follow the full GNU-standard
# we need automake 1.5
AUTOMAKE_OPTIONS = foreign 1.5

# pass most important options when "make distcheck" is used
DISTCHECK_CONFIGURE_FLAGS = --with-dune-grid=$(DUNE_GRID_ROOT) --with-dune-fem=$(DUNE_FEM_ROOT) --with-dune-common=$(DUNE_COMMON_ROOT) --with-ag-common=$(AG_COMMON_ROOT) --with-dune-istl=$(DUNE_ISTL_ROOT) --with-dune-localfunctions=$(DUNE_LOCALFUNCTIONS_ROOT)  CXX="$(CXX)" CC="$(CC)"

include $(top_srcdir)/am/global-rules

