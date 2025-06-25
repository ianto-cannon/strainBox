#!/bin/bash
#SBATCH -J spectrum
#SBATCH -t 2-00:00:00
#SBATCH -N 1
##SBATCH --mem=30G
##module load impi
module load impi
mpirun -np 1 ./spectrum
