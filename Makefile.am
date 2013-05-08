# $Id: duneproject 5842 2010-01-20 18:48:34Z joe $

# we need the module file to be able to build via dunecontrol
EXTRA_DIST=dune.module

SUBDIRS = dune src m4  doc

if BUILD_DOCS
# TODO: set up documentation tree automatically
#SUBDIRS += doc
endif

# don't follow the full GNU-standard
# we need automake 1.5
AUTOMAKE_OPTIONS = foreign 1.5

# pass most important options when "make distcheck" is used
DISTCHECK_CONFIGURE_FLAGS = --with-dune-grid=$(DUNE_GRID_ROOT) --with-dune-fem=$(DUNE_FEM_ROOT) --with-dune-common=$(DUNE_COMMON_ROOT) --with-ag-common=$(AG_COMMON_ROOT) --with-dune-istl=$(DUNE_ISTL_ROOT) --with-dune-solvers=$(DUNE_SOLVERS_ROOT) --with-dune-pdelab=$(DUNE_PDELAB_ROOT) --with-dune-localfunctions=$(DUNE_LOCALFUNCTIONS_ROOT)  CXX="$(CXX)" CC="$(CC)"

AM_CPPFLAGS = $(AM_CPPFLAGS) \
	$(DUNEMPICPPFLAGS) \
	$(UG_CPPFLAGS) \
	$(AMIRAMESH_CPPFLAGS) \
	$(ALBERTA_CPPFLAGS) $(BOOST_CPPFLAGS) $(DUNE_CPPFLAGS) \
	$(ALUGRID_CPPFLAGS)  \
	-DGRIDDIM=$(GRIDDIM) -DGRIDTYPE=$(GRIDTYPE) -DWORLDDIM=$(GRIDDIM) -DPROBLEM_NAME=$(PROBLEM)
	
#GRIDTYPE=ALBERTAGRID
GRIDTYPE=ALUGRID_SIMPLEX
GRIDDIM=2
PROBLEM_NAMES= Toy One Two Three Four Five Six Seven Eight Nine Ten
PROBLEM=Nine

include $(top_srcdir)/am/top-rules
include $(top_srcdir)/am/global-rules

