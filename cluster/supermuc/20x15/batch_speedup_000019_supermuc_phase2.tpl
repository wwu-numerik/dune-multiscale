#!/bin/bash
##
## optional: energy policy tags
#@ energy_policy_tag = et_20x15_19N_1C
#@ minimize_time_to_solution = yes
# DO NOT USE environment = COPY_ALL
#@ job_type = MPICH
#@ class = test
#@ node = 19
#@ island_count=1
# other example
#@ tasks_per_node = 28
#@ wall_clock_limit = 0:1:30
##                    1 h 20 min 30 secs
#@ job_name = msfem_20x15_19N_1C
#@ network.MPI = sn_all,not_shared,us
#@ initialdir = $(home)/multiscale-build-phase2/dune-multiscale
#@ output = job_msfem_20x15_19N_1C_$(jobid).out
#@ error = job_msfem_20x15_19N_1C_$(jobid).err
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

NODES=19
PROCS=512
OPT="$HOME/dune-multiscale-super/dune-multiscale/parameter_files/supermuc_test \
-global.datadir $HOME/multiscale-build-phase2/dune-multiscale/speedup_n532_20x15_T1 \
-grids.macro_cells_per_dim 20 -grids.micro_cells_per_macrocell_dim 15 -threading.max_count 1 "
mpiexec ${MPI} -n 512 $BIN ${OPT}