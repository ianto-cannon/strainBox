#!/bin/bash
#SBATCH -J strainBoxWe02
#SBATCH -t 2-00:00:00
#SBATCH -n 16
module load impi
mpirun -np 16 ./dropStrain
