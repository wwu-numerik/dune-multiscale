#!/usr/bin/python

import jinja2
import sys

problems = ['Easy', 'Toy', 'One', 'Two', 'Three', 'Four', 'Five', 'Six', 'Seven', 'Eight', 'Nine', 'Ten' ]
algos = ['hmm', 'msfem']

tpl = jinja2.Template(u'''
GRIDTYPE=YASPGRID
GRIDDIM=2

# tests where program to build and program to run are equal
NORMALTESTS = {% for algo in algos %}{% for problem in problems %}test_{{algo}}_{{ problem|lower }} {% endfor %}{% endfor %}

# list of tests to run
TESTS = $(NORMALTESTS)

# programs just to build when "make check" is used
check_PROGRAMS = $(NORMALTESTS)


LDFLAGS = -lboost_filesystem -lboost_system
CXXFLAGS = $(BOOST_CPPFLAGS) $(DUNE_CPPFLAGS) \\
    -DGRIDDIM=$(GRIDDIM) -D$(GRIDTYPE)

{% for algo in algos %}{% for problem in problems %}
test_{{algo}}_{{ problem|lower }}_SOURCES = $(top_srcdir)/src/elliptic_{{  algo|lower }}.cc
test_{{algo}}_{{ problem|lower }}_CPPFLAGS = -DPROBLEM_NAME={{ problem }}.cc
{% endfor %}{% endfor %}    

## distribution tarball
SOURCES = {% for algo in algos %}$(top_srcdir)/src/elliptic_{{  algo|lower }}.cc {% endfor %}

# gridcheck not used explicitly, we should still ship it :)
EXTRA_DIST = $(SOURCES)

CLEANFILES = *.gcda *.gcno semantic.cache simplex-testgrid*.dgf.* cube-testgrid*.dgf.* dgfparser.log

include $(top_srcdir)/am/global-rules

''')
fn = sys.argv[0].replace('.py','')
with open(fn, 'wb') as out:
    out.write(tpl.render(problems=problems, algos=algos))
        
