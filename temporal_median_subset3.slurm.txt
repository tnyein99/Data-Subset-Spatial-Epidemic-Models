#! /bin/bash
#-----------------------------------------
#SBATCH --job-name=temporal_median_subset3
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=10gb
#SBATCH --time=5-0:00:00
#-----------------------------------------
module load R/4.2.0

Rscript temporal_median_subset3.R  
