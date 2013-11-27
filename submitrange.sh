#!/bin/bash

for (( c=$1; c<=$2; c*=2 ))
do
    numnodes=$((($c+12-1)/12))
    qsub -l walltime=00:00:30,nodes=$numnodes:westmere:ppn=12 -t $c ../test_scaleup.cmd

done
