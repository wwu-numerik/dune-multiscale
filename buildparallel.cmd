#PBS -l nice=-1,walltime=02:26:00,nodes=11:ppn=12
#PBS -A o0num
#PBS -M r_milk01@wwu.de
#PBS -m ae
#PBS -q math
#PBS -N msfem
#PBS -j oe
#PBS -o output_build.txt

# #PBS -t 32
#PBS -o data/${PBS_JOBNAME%-*}/$PBS_JOBID/logdata/runlog.log
cd $PBS_O_WORKDIR
mkdir -p data/${PBS_JOBNAME%-*}/$PBS_JOBID/logdata/

MPIHOME=/Applic.PALMA/mpi/openmpi/gcc4.7.2/1.6.3 likwid-mpirun  -d -mpi=openmpi -np $PBS_ARRAYID  -nperdomain N:12 -- "./elliptic_msfem ../parameter_files/msfem.param -global.datadir data/${PBS_JOBNAME%-*}/$PBS_JOBID"
