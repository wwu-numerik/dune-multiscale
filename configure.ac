# -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
AC_PREREQ(2.50)
DUNE_AC_INIT # gets module version from dune.module file
AM_INIT_AUTOMAKE([subdir-objects])
AC_CONFIG_SRCDIR([src/elliptic_hmm.cc])
AM_CONFIG_HEADER([config.h])


# we need no more than the standard DE-stuff
# this module depends on dune-grid dune-fem dune-common ag-common dune-istl dune-solvers dune-pdelab dune-localfunctions
# this implies checking for [dune-common], [dune-grid], [dune-istl], [dune-localfunctions], [dune-fem]
DUNE_CHECK_ALL

# implicitly set the Dune-flags everywhere
AC_SUBST(AM_CPPFLAGS, '$(DUNE_CPPFLAGS) $(DUNEMPICPPFLAGS) $(UG_CPPFLAGS) $(AMIRAMESH_CPPFLAGS) $(ALBERTA_CPPFLAGS) $(BOOST_CPPFLAGS) $(ALUGRID_CPPFLAGS)') 
AC_SUBST(AM_LDFLAGS,  '$(DUNEMPILDFLAGS) $(UG_LDFLAGS) $(AMIRAMESH_LDFLAGS) $(BOOST_LDFLAGS) $(ALBERTA_LDFLAGS) $(ALUGRID_LDFLAGS) $(DUNE_LDFLAGS)')
AC_SUBST(AM_LDADD, '$(DUNE_LIBS) $(DUNE_LIBS) $(ALUGRID_LIBS) $(ALBERTA_LIBS) $(AMIRAMESH_LIBS) $(UG_LIBS) $(BOOST_FILESYSTEM_LIB) $(BOOST_SYSTEM_LIB) $(BOOST_TIMER_LIB) $(BOOST_CHRONO_LIB) $(DUNEMPILIBS) $(DUNE_LDADD)')


AC_CONFIG_FILES([
  Makefile
  src/Makefile
  doc/Makefile
  doc/doxygen/Makefile
  doc/doxygen/Doxyfile
  dune/Makefile
  dune/multiscale/Makefile
  dune/multiscale/test/Makefile
  m4/Makefile
  dune-multiscale.pc
])
AC_OUTPUT
# finally print the summary information
DUNE_SUMMARY_ALL
