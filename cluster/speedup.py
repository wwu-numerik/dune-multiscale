#!/usr/bin/env python

"""Simple demonstration of solving the Poisson equation in 2D using pyMOR's builtin discretizations.

Usage:
    speedup.py TEMPLATE STARTNODE POWER MACRO MICRO

Arguments:
    TEMPLATE        Which template file to use

    STARTNODE       Number of nodes to begin speedup run with

    POWER           How often to double the number of nodes

    MACRO           How many macro cells per dimension

    MICRO           How many micro cells per dimension

Options:
    -h, --help   Show this message.
"""

from jinja2 import FileSystemLoader, Environment
from os.path import dirname, realpath
from docopt import docopt
import math

try:
    loader = FileSystemLoader(dirname(realpath(__file__)), followlinks=True)
except TypeError as t:
    loader = FileSystemLoader(dirname(realpath(__file__)))

inargs = docopt(__doc__)
tpl = Environment(loader=loader).get_template(inargs['TEMPLATE'])

args = {'THREADS': 1, 'NODES': 2, 'MACRO': 8, 'MICRO': 4, 'POWER': 3, 'STARTNODE': 4 }
for key, value in args.items():
    try:
        args[key] = int(inargs.get(key, value))
    except ValueError:
        continue

nodes = args['STARTNODE']
for i, n in enumerate([int(math.pow(nodes, i)) for i in range(args['POWER'])]):
    args['NODES'] = n
    with open('speedup_{}'.format(i), 'wb') as out:
        out.write(tpl.render(**args))

