#!/bin/sh

#SBATCH --job-name=alphafold_IDR_puller
#SBATCH --account=fc_moilab
#SBATCH --partition=savio4_htc
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=5:00:00
#SBATCH --mem=9gb
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kinsey.long@berkeley.edu

source /global/scratch/projects/fc_moilab/kinseylong/CAM5/CAM5/bin/activate
python /global/scratch/projects/fc_moilab/kinseylong/CAM5/ArabidopsisIDRpuller.py

