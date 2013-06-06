GRIDTYPE=ALBERTAGRID
GRIDDIM=2
PROBLEM_NAMES= Toy One Two Three Four Five Six Seven Eight Nine Ten Eleven
PROBLEM=Nine
	
CPPFLAGS = $(AM_CPPFLAGS) -DGRIDDIM=$(GRIDDIM) -D$(GRIDTYPE) -DWORLDDIM=$(GRIDDIM) -DPROBLEM_NAME=$(PROBLEM)

# The libraries have to be given in reverse order (most basic libraries
# last).  Also, due to some misunderstanding, a lot of libraries include the
# -L option in LDFLAGS instead of LIBS -- so we have to include the LDFLAGS
# here as well.
LDFLAGS = $(AM_LDFLAGS)
LDADD = $(AM_LDADD)
