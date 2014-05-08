#PBS -l walltime=00:06:00,nodes=1:westmere:ppn=12
#PBS -A o0num
#PBS -M s_kaul01@uni-muenster.de
#PBS -m ae
#PBS -q default
#PBS -N build_parallel
#PBS -j oe
#PBS -o output_build.txt

cd $PBS_O_WORKDIR
# source ../../local/src/tbb/build/linux_intel64_gcc_cc4.7.2_libc2.5_kernel2.6.18_debug/tbbvars.sh
# export LD_LIBRARY_PATH=/Applic.PALMA/compiler/gcc/4.7.2/lib64:$LD_LIBRARY_PATH
# source ../../BOOST_PATH.sh
make -j 12 elliptic_msfem
/scratch/tmp/s_kaul01/external/bin/hpcstruct ./elliptic_msfem