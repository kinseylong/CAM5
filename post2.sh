#!/bin/sh

#SBATCH --job-name=IDR_distances
#SBATCH --account=co_moilab
#SBATCH --partition=savio3_htc
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=10:00:00
#SBATCH --mem=9gb

source /global/scratch/projects/fc_moilab/kinseylong/CAM5/CAM5/bin/activate
python /global/scratch/projects/fc_moilab/kinseylong/CAM5/IDR_distance_calc.py