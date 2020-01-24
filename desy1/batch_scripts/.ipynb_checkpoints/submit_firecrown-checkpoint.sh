#!/bin/bash   
#SBATCH -N 2
#SBATCH --qos=premium
#SBATCH --time=0:20:00
#SBATCH --job-name=FireCrown_Test
#SBATCH --image=docker:joezuntz/txpipe_cosmosis_firecrown:latest
#SBATCH --license=SCRATCH
#SBATCH --constraint=knl
#SBATCH --mail-user=elp25@duke.edu

export OMP_NUM_THREADS=1

source setup-firecrown


srun -n  64 shifter --volume=$PWD/cosmosis:/opt/cosmosis --volume $PWD/FireCrown/:/opt/firecrown --volume $PWD/TXPipe:/opt/txpipe firecrown compute /global/cscratch1/sd/elp25/FireCrown/examples/cosmicshear/cosmicshear.yaml 