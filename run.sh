#!/bin/bash
#SBATCH -J strainBox
#SBATCH -t 2-00:00:00
#SBATCH -n 16
module load impi
mpirun -np 1 ./strainBox
