
SUBDIRS =

noinst_PROGRAMS = elliptic_hmm elliptic_msfem

GRIDTYPE=ALBERTAGRID
GRIDDIM=2
HMM_PROBLEM=Easy
MSFEM_PROBLEM=Ten

elliptic_hmm_SOURCES = elliptic_hmm.cc

elliptic_msfem_SOURCES = elliptic_msfem.cc

elliptic_hmm_CPPFLAGS = $(ALL_PKG_CPPFLAGS) $(AM_CPPFLAGS) \
	$(DUNEMPICPPFLAGS) \
	$(UG_CPPFLAGS) \
	$(AMIRAMESH_CPPFLAGS) \
	$(ALBERTA_CPPFLAGS) \
	$(ALUGRID_CPPFLAGS) -DPROBLEM_NAME=$(HMM_PROBLEM) \
	-DGRIDDIM=$(GRIDDIM) -DGRIDTYPE=$(GRIDTYPE) -DWORLDDIM=$(GRIDDIM)

elliptic_msfem_CPPFLAGS = $(ALL_PKG_CPPFLAGS) $(AM_CPPFLAGS) \
	$(DUNEMPICPPFLAGS) \
	$(UG_CPPFLAGS) \
	$(AMIRAMESH_CPPFLAGS) \
	$(ALBERTA_CPPFLAGS) \
	$(ALUGRID_CPPFLAGS) -DPROBLEM_NAME=$(MSFEM_PROBLEM) \
	-DGRIDDIM=$(GRIDDIM) -D$(GRIDTYPE)

# The libraries have to be given in reverse order (most basic libraries
# last).  Also, due to some misunderstanding, a lot of libraries include the
# -L option in LDFLAGS instead of LIBS -- so we have to include the LDFLAGS
# here as well.
elliptic_hmm_LDADD = \
	$(DUNE_LDFLAGS) $(DUNE_LIBS) \
	$(ALUGRID_LDFLAGS) $(ALUGRID_LIBS) \
	$(ALBERTA_LDFLAGS) $(ALBERTA_LIBS) \
	$(AMIRAMESH_LDFLAGS) $(AMIRAMESH_LIBS) \
	$(UG_LDFLAGS) $(UG_LIBS) -lboost_filesystem -lboost_system \
	$(DUNEMPILIBS)	\
	$(LDADD)
elliptic_hmm_LDFLAGS = $(AM_LDFLAGS) \
	$(DUNEMPILDFLAGS) \
	$(UG_LDFLAGS) \
	$(AMIRAMESH_LDFLAGS) \
	$(ALBERTA_LDFLAGS) \
	$(ALUGRID_LDFLAGS) \
	$(DUNE_LDFLAGS)

elliptic_msfem_LDADD = \
	$(DUNE_LDFLAGS) $(DUNE_LIBS) \
	$(ALUGRID_LDFLAGS) $(ALUGRID_LIBS) \
	$(ALBERTA_LDFLAGS) $(ALBERTA_LIBS) \
	$(AMIRAMESH_LDFLAGS) $(AMIRAMESH_LIBS) \
	$(UG_LDFLAGS) $(UG_LIBS) -lboost_filesystem -lboost_system \
	$(DUNEMPILIBS)	\
	$(LDADD)
elliptic_msfem_LDFLAGS = $(AM_LDFLAGS) \
	$(DUNEMPILDFLAGS) \
	$(UG_LDFLAGS) \
	$(AMIRAMESH_LDFLAGS) \
	$(ALBERTA_LDFLAGS) \
	$(ALUGRID_LDFLAGS) \
	$(DUNE_LDFLAGS)

# don't follow the full GNU-standard
# we need automake 1.5
AUTOMAKE_OPTIONS = foreign 1.5

# pass most important options when "make distcheck" is used
DISTCHECK_CONFIGURE_FLAGS = --with-dune-grid=$(DUNE_GRID_ROOT) --with-dune-fem=$(DUNE_FEM_ROOT) --with-dune-common=$(DUNE_COMMON_ROOT) --with-ag-common=$(AG_COMMON_ROOT) --with-dune-istl=$(DUNE_ISTL_ROOT) --with-dune-solvers=$(DUNE_SOLVERS_ROOT) --with-dune-pdelab=$(DUNE_PDELAB_ROOT) --with-dune-localfunctions=$(DUNE_LOCALFUNCTIONS_ROOT)  CXX="$(CXX)" CC="$(CC)"

include $(top_srcdir)/am/global-rules

