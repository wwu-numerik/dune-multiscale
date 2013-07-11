GRIDTYPE=SPGRID
GRIDDIM=2
	
CPPFLAGS = $(AM_CPPFLAGS) -DGRIDDIM=$(GRIDDIM) -D$(GRIDTYPE) -DWORLDDIM=$(GRIDDIM)

# The libraries have to be given in reverse order (most basic libraries
# last).  Also, due to some misunderstanding, a lot of libraries include the
# -L option in LDFLAGS instead of LIBS -- so we have to include the LDFLAGS
# here as well.
LDFLAGS = $(AM_LDFLAGS)
LDADD = $(AM_LDADD)
