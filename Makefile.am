include $(top_srcdir)/common.mk

# we need the module file to be able to build via dunecontrol
EXTRA_DIST=dune.module

SUBDIRS = . src m4 doc dune

if BUILD_DOCS
# TODO: set up documentation tree automatically
#SUBDIRS += doc
endif

# don't follow the full GNU-standard
# we need automake 1.5
AUTOMAKE_OPTIONS = foreign 1.5

# pass most important options when "make distcheck" is used
DISTCHECK_CONFIGURE_FLAGS = --with-dune-grid=$(DUNE_GRID_ROOT) --with-dune-fem=$(DUNE_FEM_ROOT) --with-dune-common=$(DUNE_COMMON_ROOT) --with-ag-common=$(AG_COMMON_ROOT) --with-dune-istl=$(DUNE_ISTL_ROOT) --with-dune-solvers=$(DUNE_SOLVERS_ROOT) --with-dune-pdelab=$(DUNE_PDELAB_ROOT) --with-dune-localfunctions=$(DUNE_LOCALFUNCTIONS_ROOT)  CXX="$(CXX)" CC="$(CC)"




include $(top_srcdir)/am/top-rules
include $(top_srcdir)/am/global-rules

