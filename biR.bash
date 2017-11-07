#!/bin/bash
n128=0
dir0='biRelation'
theme='test'
dir=$dir0'/'$theme
if [ -d "./$dir" ]; then
    rm -r $dir/*
else
    mkdir $dir
fi
if [ "$n128" == "0" ]; then
    cp base/n128.py_pas $dir/n128.py
else
    cp base/n128.py_active $dir/n128.py
fi
cp $dir0/confirmK.py $dir
cp $dir0/confirmK.slurm $dir
cp $dir0/drawCK.py $dir
cp base/neuroAlter.py $dir
cp -R x86_64 $dir/
cd $dir
sbatch confirmK.slurm
