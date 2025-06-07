#!/bin/bash
#SBATCH -J strainBox
#SBATCH -t 2-00:00:00
#SBATCH -N 1
##SBATCH --cpus-per-task=1
module load impi
./strainBox
