# grid type, dim, etc are all defined in common.mk!
# include not needed for CMake
# include $(top_srcdir)/common.mk
###############################################################################
##########################            Libs       ##############################
###############################################################################
set(lib_LTLIBRARIES  libmultiscale_problems.la libmultiscale_common.la libmultiscale_msfem.la libmultiscale_hmm.la libmultiscale_fem.la libmultiscale_lodm.la)
######################################
##########  Lib: Common  #############
######################################
add_executable(libmultiscale_common_la ../dune/multiscale/tools/misc/outputparameter.cc
        ../dune/multiscale/hmm/cell_problem_numbering.cc
        ../dune/multiscale/hmm/cell_problem_solver.cc
        ../dune/multiscale/hmm/discrete_cell_operator.cc
        ../dune/multiscale/msfem/fem_solver.cc
        ../dune/multiscale/msfem/elliptic_msfem_matrix_assembler.cc
        ../dune/multiscale/msfem/localproblems/subgrid-list.cc
        ../dune/multiscale/msfem/msfem_grid_specifier.cc
        ../dune/multiscale/msfem/localproblems/localoperator.cc
        ../dune/multiscale/msfem/localproblems/localproblemsolver.cc
        ../dune/multiscale/common/error_calc.cc
        ../dune/multiscale/msfem/localproblems/weighted-clement-operator.cc
        ../dune/multiscale/msfem/localproblems/localsolutionmanager.cc)
######################################
##########  Lib: Problems  #############
######################################
add_executable(libmultiscale_problems_la ../dune/multiscale/problems/elliptic/eight.cc
        ../dune/multiscale/problems/elliptic/eleven.cc
        ../dune/multiscale/problems/elliptic/example.cc
        ../dune/multiscale/problems/elliptic/five.cc
        ../dune/multiscale/problems/elliptic/four.cc
        ../dune/multiscale/problems/elliptic/nine.cc
        ../dune/multiscale/problems/elliptic/one.cc
        ../dune/multiscale/problems/elliptic/seven.cc
        ../dune/multiscale/problems/elliptic/six.cc
        ../dune/multiscale/problems/elliptic/ten.cc
        ../dune/multiscale/problems/elliptic/thirteen.cc
        ../dune/multiscale/problems/elliptic/three.cc
        ../dune/multiscale/problems/elliptic/toy.cc
        ../dune/multiscale/problems/elliptic/twelve.cc
        ../dune/multiscale/problems/elliptic/two.cc
        ../dune/multiscale/problems/elliptic/spe10.cc
        ../dune/multiscale/problems/selector.cc
        ../dune/multiscale/problems/base.cc)
######################################
##########  Lib: MsFEM  ##############
######################################
add_executable(libmultiscale_msfem_la ../dune/multiscale/msfem/algorithm.cc
        ../dune/multiscale/msfem/msfem_solver.cc)
######################################
###########  Lib: HMM  ###############
######################################
add_executable(libmultiscale_hmm_la ../dune/multiscale/hmm/algorithm.cc
        ../dune/multiscale/hmm/algorithm_error.cc
        ../dune/multiscale/hmm/algorithm_step.cc
        ../dune/multiscale/hmm/elliptic_hmm_matrix_assembler.cc)
######################################
###########  Lib: FEM  ###############
######################################
add_executable(libmultiscale_fem_la ../dune/multiscale/fem/algorithm.cc)
######################################
###########  Lib: LODM  ###############
######################################
add_executable(libmultiscale_lodm_la ../dune/multiscale/msfem/rigorous.cc
        ../dune/multiscale/msfem/rigorous_msfem_solver.cc)
###############################################################################
##########################        Executables    ##############################
###############################################################################
set(noinst_PROGRAMS  elliptic_hmm elliptic_msfem elliptic_fem elliptic_hom_fem elliptic_lodm)
set(SOURCE_PREFIX elliptic_hmm elliptic_msfem elliptic_fem elliptic_hom_fem elliptic_lodm)
######################################
##############  HMM  #################
######################################
add_executable(elliptic_hmm elliptic_hmm.cc)
set(elliptic_hmm_DEPENDENCIES  libmultiscale_problems.la libmultiscale_common.la libmultiscale_hmm.la)
# The following lines are comment out as the translation is not clear. TODO review
#
# set(elliptic_hmm_ LDFLAGS${LDFLAGS})
# The following lines are comment out as the translation is not clear. TODO review
#
# set(elliptic_hmm_ LDADD-lmultiscale_hmm -lmultiscale_problems -lmultiscale_common ${LDADD})
######################################
##############  MSFEM  ###############
######################################
add_executable(elliptic_msfem elliptic_msfem.cc)
set(elliptic_msfem_DEPENDENCIES  libmultiscale_problems.la libmultiscale_common.la libmultiscale_msfem.la)
set(elliptic_msfem_LDPFLAGS  ${LDPFLAGS})
# The following lines are comment out as the translation is not clear. TODO review
#
# set(elliptic_msfem_ LDADD-lmultiscale_msfem -lmultiscale_problems -lmultiscale_common ${LDADD})
######################################
###########  LOD Method  #############
######################################
add_executable(elliptic_lodm elliptic_lodm.cc)
set(elliptic_lodm_DEPENDENCIES  libmultiscale_problems.la libmultiscale_common.la libmultiscale_lodm.la)
set(elliptic_lodm_LDPFLAGS  ${LDPFLAGS})
# The following lines are comment out as the translation is not clear. TODO review
#
# set(elliptic_lodm_ LDADD-lmultiscale_lodm -lmultiscale_problems -lmultiscale_common ${LDADD})
######################################
##############  FEM  #################
######################################
add_executable(elliptic_fem elliptic_fem.cc)
set(elliptic_fem_DEPENDENCIES  libmultiscale_problems.la libmultiscale_common.la libmultiscale_fem.la)
# The following lines are comment out as the translation is not clear. TODO review
#
# set(elliptic_fem_ LDFLAGS${LDFLAGS})
# The following lines are comment out as the translation is not clear. TODO review
#
# set(elliptic_fem_ LDADD-lmultiscale_fem -lmultiscale_problems -lmultiscale_common ${LDADD})
######################################
############  HOM FEM  ###############
######################################
add_executable(elliptic_hom_fem elliptic_hom_fem.cc)
set(elliptic_hom_fem_DEPENDENCIES  libmultiscale_problems.la libmultiscale_common.la)
# The following lines are comment out as the translation is not clear. TODO review
#
# set(elliptic_hom_fem_ LDFLAGS${LDFLAGS})
# The following lines are comment out as the translation is not clear. TODO review
#
# set(elliptic_hom_fem_ LDADD-lmultiscale_fem -lmultiscale_problems -lmultiscale_common ${LDADD})
###############################################################################
##########################            Misc       ##############################
###############################################################################
# don't follow the full GNU-standard
# we need automake 1.5
set(AUTOMAKE_OPTIONS  foreign 1.5)
# pass most important options when "make distcheck" is used
set(DISTCHECK_CONFIGURE_FLAGS  --with-dune-grid=${DUNE_GRID_ROOT} --with-dune-fem=${DUNE_FEM_ROOT} --with-dune-common=${DUNE_COMMON_ROOT} --with-ag-common=${AG_COMMON_ROOT} --with-dune-istl=${DUNE_ISTL_ROOT} --with-dune-localfunctions=${DUNE_LOCALFUNCTIONS_ROOT}  CXX="${CXX}" CC="${CC}")
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
