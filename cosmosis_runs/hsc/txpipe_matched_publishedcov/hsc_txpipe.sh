#!/bin/bash -l
#SBATCH -e slurm.err
#SBATCH -N 1
#SBATCH -n 40
#SBATCH --time 10:00:00
#SBATCH --mail-user=elp25@duke.edu
#SBATCH --mail-type=ALL
#SBATCH --mem-per-cpu=4G
#SBATCH -p cosmology
bash
#source /hpc/group/cosmology/loadEnv.sh
source /hpc/group/cosmology/cosmosisEnv.sh /hpc/group/cosmology/cosmosis_env
export OMP_NUM_THREADS=1
export SLURM_WHOLE=1
#export OMP_PLACES=threads
#export OMP_PROC_BIND=spread

mpirun -n $SLURM_NTASKS cosmosis --mpi hsc_params.ini
