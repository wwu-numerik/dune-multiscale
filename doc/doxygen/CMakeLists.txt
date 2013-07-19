
set(BASEDIR ../..)
set(CURDIR doc/doxygen)
# include not needed for CMake
# include $(top_srcdir)/am/doxygen
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
# Create Doxyfile.in and Doxyfile, and doxygen documentation
add_doxygen_target()