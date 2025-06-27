#!/bin/bash
#SBATCH -J strainBoxWe05
#SBATCH -t 2-00:00:00
#SBATCH -n 32
module load impi
mpirun -np 32 ./dropStrain
