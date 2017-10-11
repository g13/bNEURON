#!/bin/bash
n128=0
theme='mpiTest-active-debug'

if [ -d "./$theme" ]; then
    rm -r $theme/*
else
    mkdir $theme
fi

if [ "$n128" == "0" ]; then
    cp n128.py_pas $theme/n128.py
else
    cp n128.py_active $theme/n128.py
fi
cp getNib_mpi.slurm $theme
cp getNib_mpi_0.slurm $theme
cp neuroAlter.py $theme
cp getNib_mpi.py $theme
cp -R x86_64 $theme/

cd $theme
python -m py_compile getNib_mpi.py

sbatch --export=theme=$theme getNib_mpi.slurm
sbatch --export=theme=$theme getNib_mpi_0.slurm
