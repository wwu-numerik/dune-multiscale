#!/usr/bin/env python

tpl = '''#!/bin/bash -l 
#SBATCH --nodes={{ NODES }}
#SBATCH --ntasks-per-node={{ PROCS }}
#SBATCH --ntasks={{ NODES * PROCS }}
#SBATCH --cpus-per-task=1
#SBATCH --mem=4gb 
#SBATCH --time=02:00:00 
#SBATCH --account=UniKoeln 
#SBATCH --partition=mpi
#SBATCH --constraint=[inca12*{{ NODES }}]
#SBATCH --exclusive
#SBATCH -o out_{{ NODES }}_{{ PROCS }}_{{ THREADS }}_{{ MACRO }}x{{ MICRO }}_%j_%N.log
#SBATCH -e out_{{ NODES }}_{{ PROCS }}_{{ THREADS }}_{{ MACRO }}x{{ MICRO }}_%j_%N.err
#SBATCH --mail-user=rmilk@uni-koeln.de
#SBATCH --mail-type=end 

source $HOME/.bashrc
source $HOME/.modules
MA={{ MACRO }}
MI={{ MICRO }}
BIN=/scratch/rmilk/multiscale-build/dune-multiscale/elliptic_msfem
COMMON="-threading.max_count {{ THREADS }} -threading.partition_factor 1"
DATADIR=speedup
SUF="_T01"

OPT="/home/rmilk/cmake_multiscale/dune-multiscale/parameter_files/cheops.ini -global.datadir /scratch/rmilk/dune/single_run_np{{ NODES }} -grids.macro_cells_per_dim {{ MACRO }} -grids.micro_cells_per_macrocell_dim {{ MICRO }}  ${COMMON}"
srun -n $PROCS -N $NODES $BIN ${OPT}
'''

from jinja2 import Environment as env

args = { 'PROCS' : 4, 'NODES': 8, 'THREADS': 1, 'MACRO': 8, 'MICRO': 4}
kk = env().from_string(tpl)
print(kk.render(**args))

