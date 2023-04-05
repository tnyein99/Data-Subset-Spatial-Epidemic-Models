#! /bin/bash
#-----------------------------------------
#SBATCH --job-name=k_per_corner
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=5gb
#SBATCH --time=5-0:00:00
#-----------------------------------------
module load R/4.2.0

Rscript k_per_corner.R  




