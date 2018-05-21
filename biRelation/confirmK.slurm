#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10G
#SBATCH --time=01:00:00
#SBATCH -a 1-4
#SBATCH -o cK_%A_%a.out
#SBATCH --job-name=cK
#SBATCH --mail-type=END,FAIL

module purge
module load anaconda
module list

date
source activate py2
v0=-74
dV=4
seed=238660
dt=0
fmt='eps'
#dtRange = np.array([0,4,8,12,20,30,70,130,210],dtype='float')
python confirmK.py --dV=$dV --vmin=$v0 --vid=$(($SLURM_ARRAY_TASK_ID-1)) --seed=$seed --dt=$dt --fmt=$fmt --so=1 --nf=1

date
hostname
