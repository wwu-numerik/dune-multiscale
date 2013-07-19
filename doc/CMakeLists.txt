
set(SUBDIRS  doxygen)
set(CURDIR  doc)
set(BASEDIR  ..)
# add list of html files to generate from wml
set(PAGES )
set(docdir ${CMAKE_INSTALL_DOCDIR}/dune-multiscale)
# include not needed for CMake
# include $(top_srcdir)/am/webstuff
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
foreach(i ${SUBDIRS})
  if(${i} STREQUAL "test")
    set(opt EXCLUDE_FROM_ALL)
  endif(${i} STREQUAL "test")
  add_subdirectory(${i} ${opt})
endforeach(i ${SUBDIRS})
