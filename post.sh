#!/bin/sh

#SBATCH --job-name=alphafold_IDR_puller

#SBATCH --account=fc_moilab

#SBATCH --partition=savio3_htc

#SBATCH --ntasks-per-node=1

#SBATCH --cores=4

#SBATCH --time=24:00:00 --mem 9gb

source /global/scratch/projects/fc_moilab/kinseylong/CAM5/CAM5/bin/activate
python /global/scratch/projects/fc_moilab/kinseylong/CAM5/ArabidopsisIDRpuller.py

