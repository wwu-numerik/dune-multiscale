#!/bin/bash
##
## optional: energy policy tags
#@ energy_policy_tag = et_{{ MACRO }}x{{ MICRO }}_{{ NODES }}N_{{ THREADS }}C
#@ minimize_time_to_solution = yes
# DO NOT USE environment = COPY_ALL
#@ job_type = MPICH
#@ class = {{ 'test' if NODES < 20 else 'general' }}
#@ node = {{ NODES }}
#@ island_count=1
# other example
#@ tasks_per_node = 28
#@ wall_clock_limit = 0:2:30
##                    1 h 20 min 30 secs
#@ job_name = msfem_{{ MACRO }}x{{ MICRO }}_{{ NODES }}N_{{ THREADS }}C
#@ network.MPI = sn_all,not_shared,us
#@ initialdir = $(home)/multiscale-build-phase2/dune-multiscale
#@ output = job_msfem_{{ MACRO }}x{{ MICRO }}_{{ NODES }}N_{{ THREADS }}C_$(jobid).out
#@ error = job_msfem_{{ MACRO }}x{{ MICRO }}_{{ NODES }}N_{{ THREADS }}C_$(jobid).err
#@ notification=always
#@ notify_user=rene.milk@wwu.de
#@ queue
. /etc/profile
. /etc/profile.d/modules.sh
#setup of environment
source $HOME/.modules
#optional: 
#module load mpi_pinning/hybrid_blocked

BIN=$HOME/multiscale-build-phase2/dune-multiscale/elliptic_msfem
# MPI="-check-mpi -prepend-rank"
# MPI="-prepend-rank"
# MPI="-binding pin=1;cell=unit;map=spread"

NODES={{ NODES }}
PROCS={{ POWER2_28 }}
OPT="$HOME/dune-multiscale-super/dune-multiscale/parameter_files/supermuc_test \
-global.datadir $HOME/multiscale-build-phase2/dune-multiscale/speedup_n{{ 28 * NODES }}_{{ MACRO }}x{{ MICRO }}_T{{ THREADS }} \
-grids.macro_cells_per_dim {{ MACRO }} -grids.micro_cells_per_macrocell_dim {{ MICRO }} -threading.max_count {{ THREADS }} "
mpiexec ${MPI} -n {{ POWER2_28 }} $BIN ${OPT}
