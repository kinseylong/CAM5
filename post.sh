#!/bin/sh

#SBATCH --job-name=alphafold_IDR_puller

#SBATCH --account=fc_moilab

#SBATCH --partition=savio3_htc

#SBATCH --ntasks-per-node=1

#SBATCH --cores=1

#SBATCH --time=48:00:00 --mem 9gb