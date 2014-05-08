#PBS -l walltime=00:00:20,nodes=2:westmere:ppn=12
#PBS -A o0num
#PBS -M s_kaul01@uni-muenster.de
#PBS -m ae
#PBS -q math
#PBS -N test_scaleup_128_466
#PBS -j oe
#PBS -t 16
#PBS -o data/${PBS_JOBNAME%-*}/$PBS_JOBID/logdata/runlog.log
cd $PBS_O_WORKDIR
mkdir -p data/${PBS_JOBNAME%-*}/$PBS_JOBID/logdata/
# export LD_LIBRARY_PATH=/Applic.PALMA/compiler/gcc/4.7.2/lib64:$LD_LIBRARY_PATH
### ${PBS_JOBNAME%-*} deletes "-" and the array-id from the jobname

ulimit -q 13107200
time mpiexec -np $PBS_ARRAYID ./elliptic_msfem ../build/msfem.param -global.datadir data/${PBS_JOBNAME%-*}/$PBS_JOBID


# mpiexec -np $PBS_ARRAYID /scratch/tmp/s_kaul01/external/bin/hpcrun -o data/${PBS_JOBNAME%-*}/$PBS_JOBID/hpctool_measurements ./elliptic_msfem ../build/msfem.param -global.datadir data/${PBS_JOBNAME%-*}/$PBS_JOBID 
# /scratch/tmp/s_kaul01/external/bin/hpcprof -M stats -S elliptic_msfem.hpcstruct -I ../dune/+ -I ../src/+ -I ../../dune-fem/dune/+ -I ../../dune-spgrid/dune/+ -I ../../dune-common/dune/+ -I ../../dune-grid/dune/+ -I ../../dune-localfunctions/dune/+ -I ../../dune-stuff/dune/+ -I ../../dune-geometry/dune/+ -I ../../dune-istl/dune/+ -o data/${PBS_JOBNAME%-*}/$PBS_JOBID/hpctool_database data/${PBS_JOBNAME%-*}/$PBS_JOBID/hpctool_measurements/