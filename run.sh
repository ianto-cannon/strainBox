#!/bin/bash
#SBATCH -J spectrum
#SBATCH -t 2-00:00:00
#SBATCH -N 1
#SBATCH --exclude=node15,node16,node01,node02,node08
##SBATCH --mem=30G
##module load impi
module load impi
mpirun -np 1 ./spectrum
