#!/bin/bash
#SBATCH -J strainBox
#SBATCH -t 2-00:00:00
#SBATCH -n 1
module load impi
mpirun -np 1 ./dropStrain
