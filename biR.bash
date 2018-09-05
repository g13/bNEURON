#!/bin/bash
n128=2
dir0='biRelation'
theme='pasT_plot'
nc=1
dir=$dir0'/'$theme
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
cp $dir0/confirmK.py $dir
cp $dir0/confirmK.slurm $dir
cp $dir0/drawCK.py $dir
cp base/neuroAlter.py $dir
cp -R x86_64 $dir/
cd $dir
export nc
sbatch --export=ALL confirmK.slurm
