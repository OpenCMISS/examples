#!/bin/bash
#SBATCH -J Monodomain
#SBATCH -A uoa00379  
#SBATCH -C avx
#SBATCH --time=00:30:00 
#SBATCH --ntasks=32
#SBATCH --mem-per-cpu=1024  
#SBATCH --workdir=/projects/uoa00379/OpenCMISS-mpch/examples/Bioelectrics/build/Monodomain

module load intel/2015.02

time srun ./monodomain.x 0.001 1 250 n98.xml 125
