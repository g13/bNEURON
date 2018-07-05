#!/bin/bash
n128=0
dir0=library
nc=1
theme='pasT_nc'
dir=$dir0/$theme

if [ -d "./$dir" ]; then
    rm -r $dir/*
else
    mkdir $dir
fi

if [ "$n128" == "0" ]; then
    cp base/n128.py_pas $dir/n128.py
fi
if [ "$n128" == "1" ]; then
    cp base/n128.py_active $dir/n128.py
fi
if [ "$n128" == "2" ]; then
    cp base/n128.py_simple $dir/n128.py
fi
cp $dir0/getNib_mpi.py $dir
cp $dir0/getNib_mpi.slurm $dir
cp $dir0/getNib_mpi_0.slurm $dir
cp base/neuroAlter.py $dir
cp -r x86_64 $dir

cd $dir
python -m py_compile getNib_mpi.py
export theme
export nc
sbatch --export=ALL getNib_mpi.slurm
sbatch --export=ALL getNib_mpi_0.slurm
