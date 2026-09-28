#!/bin/bash
#SBATCH --job-name=plotfin
#SBATCH --partition=astro_long
#SBATCH --nodes=1
#SBATCH --mem=10G
#SBATCH --mail-type=ALL
 
module load astro matlab
unset DISPLAY
cd /groups/astro/ardash/Documents/masstest/plot/matlab/
matlab -nodisplay < plotFinal.m > outfile.txt
