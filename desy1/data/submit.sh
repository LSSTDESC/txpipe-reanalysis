#!/bin/bash -l

#SBATCH -C haswell
#SBATCH -q regular
#SBATCH -N 4
#SBATCH -t 20:00:00
#SBATCH -A m1727
#SBATCH -J desy1
#SBATCH --mail-user=chihway@kicp.uchicago.edu
#SBATCH --mail-type=ALL
#SBATCH --output=cosmic_shear_desy1.out
#SBATCH --error=cosmic_shear_desy1.err

export OMP_NUM_THREADS=1

source /global/homes/c/chihway/cosmosis_y1/cosmosis/config/setup-cosmosis-nersc

mpirun -n 32 cosmosis params.ini

