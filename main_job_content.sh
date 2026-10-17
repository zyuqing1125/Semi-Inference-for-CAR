#!/bin/bash
#SBATCH --job-name=new_bootstrap
#SBATCH --time=00:30:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5G

module load R

Rscript job_run.R $1
