#!/bin/bash

mode=1
n128=1
theme='freqLimit-n'$n128-'m'$mode-'p'
## p-preset; r-random

if [ -d "./$theme" ]; then
    rm -r $theme/*
else
    mkdir $theme
fi
cp freqLimit.slurm $theme
if [ "$mode" == "0" ]; then
    cp freqLimit.py $theme
else
    cp freqLimit4.py $theme
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

sbatch --export=mode=$mode freqLimit.slurm
