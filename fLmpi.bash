#!/bin/bash

mode=0
n128=0
theme='fL10-n'$n128-'m'$mode-p
## p-preset; r-random

if [ -d "./$theme" ]; then
    rm -r $theme/*
else
    mkdir $theme
fi
cp freqMpi.slurm $theme
if [ "$mode" == "0" ]; then
    cp freqLimit0_mpi.py $theme
else
    cp freqLimit_mpi.py $theme
fi
if [ "$n128" == "0" ]; then
    cp n128.py_pas $theme/n128.py
else
    cp n128.py_active $theme/n128.py
fi
cp neuroAlter.py $theme
cp getPSP.py $theme
cp -R x86_64 $theme/

cd $theme

sbatch --export=mode=$mode freqMpi.slurm
