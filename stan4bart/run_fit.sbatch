#!/bin/bash
#
#SBATCH --job-name=s4b_sim
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --array=1-1000 # this creates an array!
#SBATCH --time=330:00
#SBATCH --mem=10GB
#SBATCH --output=stan4bart_%A_%a.out
#SBATCH --error=stan4bart_%A_%a.err
module purge
module load r/gcc/4.1.2



Rscript fit.R
