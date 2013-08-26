/* begin dune-multiscale */

/* Define to the version of dune-multiscale */
#define DUNE_MULTISCALE_VERSION "${DUNE_MULTISCALE_VERSION}"

/* Define to the major version of dune-multiscale */
#define DUNE_MULTISCALE_VERSION_MAJOR ${DUNE_MULTISCALE_VERSION_MAJOR}

/* Define to the minor version of dune-multiscale */
#define DUNE_MULTISCALE_VERSION_MINOR ${DUNE_MULTISCALE_VERSION_MINOR}

/* Define to the revision of dune-multiscale */
#define DUNE_MULTISCALE_VERSION_REVISION ${DUNE_MULTISCALE_VERSION_REVISION}

#undef COMMIT

#ifndef DUNE_MULTISCALE_CONFIG_H
#define DUNE_MULTISCALE_CONFIG_H

/* the git tag / commit we build from */
#define COMMIT "@COMMIT@"

#define ALBERTA_DIM WORLDDIM
#define ALBERTA_DEBUG 1
#cmakedefine @ENABLE_ALUGRID@ 
#cmakedefine @ENABLE_ABERTA@
#define ENABLE_MPI @ENABLE_MPI@

#define DUNE_COMMON_FIELDVECTOR_SIZE_IS_METHOD 1

#define HAVE_PETSC ENABLE_PETSC
#ifdef NDEBUG
  #define DNDEBUG
#endif

//compiler quirks workarounds
#ifdef __clang__
  class type_info;
  #define BOOST_HAS_RVALUE_REFS 1
#endif

#endif  /* DUNE_MULTISCALE_CONFIG_H */

/* end dune-multiscale */