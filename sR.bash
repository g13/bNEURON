#!/bin/bash
n128=0
theme='massive-spike'
lib='mpiTest-massive-NEURON.mat'
rE=40
rI=60
t=1000
v0=1
tref=13
seed=193864

datafdr='../../data'
dir0='singleRun'
dir1=$dir0'/'$theme
lib=$datafdr'/'$lib
if [ -d "$dir1" ]; then
    rm -r $dir1/*
else
    mkdir $dir1
fi

if [ "$n128" == "0" ]; then
    cp base/n128.py_pas $dir1/n128.py
else
    cp base/n128.py_active $dir1/n128.py
fi
cp base/neuroAlter.py $dir1
cp $dir0/sR.slurm $dir1
cp $dir0/sR_init.cfg $dir1
cp $dir0/sR $dir1
cp $dir0/bnsynCompare.m $dir1
cp -R x86_64 $dir1

cd $dir1

export lib
export v0
export theme
export rE
export rI
export t
export tref
export seed
sbatch --export=ALL sR.slurm
