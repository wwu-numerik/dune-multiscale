
# own defines
GRIDDIM=2
POLORDER=0
GRIDTYPE=ALUGRID_SIMPLEX
#GRIDTYPE=SGRID

AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) -DGRIDDIM=$(GRIDDIM) -DPOLORDER=$(POLORDER) -D$(GRIDTYPE)

noinst_PROGRAMS = disc_funcs

disc_funcs_SOURCES = disc_funcs.cc

disc_funcs_CXXFLAGS = $(MPI_CPPFLAGS)  $(AM_CPPFLAGS)
disc_funcs_LDADD = $(MPI_LDFLAGS) $(ALUGRID_LDFLAGS) $(ALUGRID_LIBS)  \
$(MPI_LDFLAGS) $(DUNE_LDFLAGS) $(DUNE_LIBS) $(GRAPE_LDFLAGS) $(GRAPE_LIBS) $(BLAS_LIBS)

# don't follow the full GNU-standard
# we need automake 1.5
AUTOMAKE_OPTIONS = foreign 1.5

# pass most important options when "make distcheck" is used
DISTCHECK_CONFIGURE_FLAGS = --with-dune=$(DUNEROOT) CXX="$(CXX)" CC="$(CC)"

include $(top_srcdir)/am/global-rules
