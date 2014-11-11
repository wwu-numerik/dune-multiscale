#PBS -l nice=-18,walltime=03:00:00,nodes=11:ppn=12
#PBS -A o0num
#PBS -k oe
#PBS -m ae
#PBS -M r_milk01@wwu.de
#PBS -q math
#PBS -N week_21_tests_S
#PBS -j oe
# #PBS -t 32
#PBS -o data/${PBS_JOBNAME%-*}/$PBS_JOBID/logdata/runlog.log
cd $PBS_O_WORKDIR
mkdir -p data/${PBS_JOBNAME%-*}/$PBS_JOBID/logdata/

MPIHOME=/Applic.PALMA/mpi/openmpi/gcc4.7.2/1.6.3 likwid-mpirun  -d -mpi=openmpi -np $PBS_ARRAYID  -nperdomain N:12 -- "./elliptic_msfem ../build/msfem.param -global.datadir data/${PBS_JOBNAME%-*}/$PBS_JOBID"

# mpiexec -np $PBS_ARRAYID ./elliptic_msfem ../build/msfem.param -global.datadir data/${PBS_JOBNAME%-*}/$PBS_JOBID
# mpirun -np $PBS_ARRAYID amplxe-cl -result-dir data/$PBS_JOBID/ampl-xe -collect hotspots -- ./elliptic_msfem ../build/msfem.param


