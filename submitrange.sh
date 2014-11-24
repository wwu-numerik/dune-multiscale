#!/bin/bash
echo "Logging-level richtig eingestellt?"

for (( c=$1; c<=$2; c*=2 ))
do
    numnodes=$((($c+12-1)/12))
    qsub -l nice=-18,walltime=03:00:00,nodes=$numnodes:westmere:ppn=12 -t $c ./speedup

done
