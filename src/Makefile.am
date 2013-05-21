# grid type, dim, etc are all defined in common.mk!
include $(top_srcdir)/common.mk


###############################################################################
##########################            Libs       ##############################
###############################################################################

lib_LTLIBRARIES = libmultiscale_problems.la libmultiscale_common.la libmultiscale_msfem.la libmultiscale_hmm.la libmultiscale_fem.la

######################################
##########  Lib: Common  #############
######################################
libmultiscale_common_la_SOURCES=../dune/multiscale/tools/misc/outputparameter.cc \
        ../dune/multiscale/hmm/cell_problem_numbering.cc \
        ../dune/multiscale/hmm/cell_problem_solver.cc \
        ../dune/multiscale/hmm/discrete_cell_operator.cc \
        ../dune/multiscale/msfem/fem_solver.cc \
        ../dune/multiscale/msfem/elliptic_msfem_matrix_assembler.cc \
        ../dune/multiscale/msfem/localproblems/subgrid-list.cc \
        ../dune/multiscale/msfem/msfem_grid_specifier.cc \
        ../dune/multiscale/msfem/localproblems/localoperator.cc \
        ../dune/multiscale/msfem/localproblems/localproblemsolver.cc
libmultiscale_common_la_CPPFLAGS=$(CPPFLAGS)

######################################
##########  Lib: Problems  #############
######################################
libmultiscale_problems_la_SOURCES=../dune/multiscale/problems/elliptic/eight.cc \
	../dune/multiscale/problems/elliptic/example.cc \
	../dune/multiscale/problems/elliptic/five.cc \
	../dune/multiscale/problems/elliptic/four.cc \
	../dune/multiscale/problems/elliptic/nine.cc \
	../dune/multiscale/problems/elliptic/one.cc \
	../dune/multiscale/problems/elliptic/seven.cc \
	../dune/multiscale/problems/elliptic/six.cc \
	../dune/multiscale/problems/elliptic/ten.cc \
	../dune/multiscale/problems/elliptic/three.cc \
	../dune/multiscale/problems/elliptic/toy.cc \
	../dune/multiscale/problems/elliptic/two.cc
libmultiscale_problems_la_CPPFLAGS=$(CPPFLAGS)

######################################
##########  Lib: MsFEM  ##############
######################################
libmultiscale_msfem_la_SOURCES=../dune/multiscale/msfem/algorithm.cc \
        ../dune/multiscale/msfem/msfem_solver.cc
libmultiscale_msfem_la_CPPFLAGS=$(CPPFLAGS)

######################################
###########  Lib: HMM  ###############
######################################
libmultiscale_hmm_la_SOURCES=../dune/multiscale/hmm/algorithm.cc \
        ../dune/multiscale/hmm/algorithm_error.cc \
        ../dune/multiscale/hmm/algorithm_step.cc \
        ../dune/multiscale/hmm/elliptic_hmm_matrix_assembler.cc
libmultiscale_hmm_la_CPPFLAGS=$(CPPFLAGS)

######################################
###########  Lib: FEM  ###############
######################################
libmultiscale_fem_la_SOURCES=../dune/multiscale/fem/algorithm.cc
libmultiscale_fem_la_CPPFLAGS=$(CPPFLAGS)

###############################################################################
##########################        Executables    ##############################
###############################################################################

noinst_PROGRAMS = elliptic_hmm elliptic_msfem elliptic_fem elliptic_hom_fem elliptic_rigorous_msfem
SOURCE_PREFIX=elliptic_hmm elliptic_msfem elliptic_fem elliptic_hom_fem elliptic_rigorous_msfem

######################################
##############  HMM  #################
######################################
elliptic_hmm_SOURCES = elliptic_hmm.cc 
elliptic_hmm_DEPENDENCIES = libmultiscale_common.la libmultiscale_hmm.la
elliptic_hmm_CPPFLAGS=$(CPPFLAGS) 
elliptic_hmm_LDFLAGS=$(LDFLAGS) 
elliptic_hmm_LDADD=-lmultiscale_hmm -lmultiscale_problems -lmultiscale_common $(LDADD)

######################################
##############  MSFEM  ###############
######################################
elliptic_msfem_SOURCES = elliptic_msfem.cc
elliptic_msfem_DEPENDENCIES = libmultiscale_common.la libmultiscale_msfem.la
elliptic_msfem_CPPFLAGS = $(CPPFLAGS)
elliptic_msfem_LDPFLAGS = $(LDPFLAGS)
elliptic_msfem_LDADD=-lmultiscale_msfem -lmultiscale_problems -lmultiscale_common $(LDADD)
        
######################################
#########  Rigorous MSFEM  ###########
######################################
elliptic_rigorous_msfem_SOURCES = elliptic_rigorous_msfem.cc \
        ../dune/multiscale/msfem/rigorous.cc \
        ../dune/multiscale/msfem/rigorous_msfem_solver.cc 
elliptic_rigorous_msfem_DEPENDENCIES = libmultiscale_common.la libmultiscale_fem.la
elliptic_rigorous_msfem_CPPFLAGS = $(CPPFLAGS)
elliptic_rigorous_msfem_LDPFLAGS = $(LDPFLAGS)
elliptic_rigorous_msfem_LDADD=-lmultiscale_fem -lmultiscale_problems -lmultiscale_common $(LDADD)

######################################
##############  FEM  #################
######################################        
elliptic_fem_SOURCES = elliptic_fem.cc 
elliptic_fem_DEPENDENCIES = libmultiscale_common.la libmultiscale_fem.la
elliptic_fem_CPPFLAGS=$(CPPFLAGS) 
elliptic_fem_LDFLAGS=$(LDFLAGS) 
elliptic_fem_LDADD=-lmultiscale_fem -lmultiscale_problems -lmultiscale_common $(LDADD)

######################################
############  HOM FEM  ###############
######################################
elliptic_hom_fem_SOURCES = elliptic_hom_fem.cc 
elliptic_hom_fem_DEPENDENCIES = libmultiscale_common.la
elliptic_hom_fem_CPPFLAGS = $(CPPFLAGS)
elliptic_hom_fem_LDFLAGS=$(LDFLAGS) 
elliptic_hom_fem_LDADD=-lmultiscale_fem -lmultiscale_problems -lmultiscale_common $(LDADD)

###############################################################################
##########################            Misc       ##############################
###############################################################################

# don't follow the full GNU-standard
# we need automake 1.5
AUTOMAKE_OPTIONS = foreign 1.5

# pass most important options when "make distcheck" is used
DISTCHECK_CONFIGURE_FLAGS = --with-dune-grid=$(DUNE_GRID_ROOT) --with-dune-fem=$(DUNE_FEM_ROOT) --with-dune-common=$(DUNE_COMMON_ROOT) --with-ag-common=$(AG_COMMON_ROOT) --with-dune-istl=$(DUNE_ISTL_ROOT) --with-dune-localfunctions=$(DUNE_LOCALFUNCTIONS_ROOT)  CXX="$(CXX)" CC="$(CC)"

include $(top_srcdir)/am/global-rules
