#!/bin/bash

n128=0
mode=0
theme='fLsp-n'$n128-'m'$mode
## p-preset; r-random

if [ -d "./$theme" ]; then
    rm -r $theme/*
else
    mkdir $theme
fi
cp fL_spList.slurm $theme
if [ "$mode" == "0" ]; then
    cp fL0_spList.py $theme
else
    cp fL_spList.py $theme
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

sbatch --export=mode=$mode fL_spList.slurm
