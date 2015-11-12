#!/usr/bin/env python

from jinja2 import Environment as env




args = { 'PROCS' : 4, 'NODES': 8, 'THREADS': 1, 'MACRO': 8, 'MICRO': 4}
kk = env().from_string(tpl)
print(kk.render(**args))

