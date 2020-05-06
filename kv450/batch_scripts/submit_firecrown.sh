#!/bin/bash -l

#SBATCH -C haswell
#SBATCH -q regular
#SBATCH -N 4
#SBATCH -t 20:00:00
#SBATCH -A m1727
#SBATCH -J firecrown_cosmicshear_kv450
#SBATCH --mail-user=elp25@duke.edu
#SBATCH --mail-type=ALL
#SBATCH --output=cosmic_shear_kv450.out
#SBATCH --error=cosmic_shear_kv450.err

export OMP_NUM_THREADS=1

source activate firecrown 

cd /global/cscratch1/sd/elp25/txpipe-cosmodc2/firecrown_config
mpirun -n 32 firecrown run-cosmosis cosmodc2_firecrown_real_fisher.yaml




