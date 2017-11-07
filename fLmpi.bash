#!/bin/bash

mode=0
n128=0
dir0='freqLimit'
theme='fL10-n'$n128-'m'$mode-p
## p-preset; r-random
dir=$dir0'/'$theme

if [ -d "./$dir" ]; then
    rm -r $dir/*
else
    mkdir $dir
fi
cp $dir0/freqMpi.slurm $dir
if [ "$mode" == "0" ]; then
    cp $dir0/freqLimit0_mpi.py $dir
else
    cp $dir0/freqLimit_mpi.py $dir
fi
if [ "$n128" == "0" ]; then
    cp base/n128.py_pas $dir/n128.py
else
    cp base/n128.py_active $dir/n128.py
fi
cp base/neuroAlter.py $dir
cp -R x86_64 $dir

cd $dir

export mode
sbatch --export=ALL freqMpi.slurm
