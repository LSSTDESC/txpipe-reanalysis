#!/bin/bash   
#SBATCH -N 2
#SBATCH --qos=regular
#SBATCH --time=15:00:00
#SBATCH --job-name=FireCrown_Test
#SBATCH --image=docker:joezuntz/txpipe_cosmosis_firecrown:latest
#SBATCH --license=SCRATCH
#SBATCH --constraint=knl
#SBATCH --mail-user=elp25@duke.edu

export OMP_NUM_THREADS=1

source setup-firecrown


srun -n  64 shifter --volume=/global/cscratch1/sd/elp25/cosmosis:/opt/cosmosis --volume /global/cscratch1/sd/elp25/FireCrown/:/opt/firecrown --volume /global/cscratch1/sd/elp25/TXPipe:/opt/txpipe firecrown compute /global/cscratch1/sd/elp25/FireCrown/examples/des_y1_3x2pt/des_y1_3x2pt.yaml 

srun -n  64 shifter --volume=/global/cscratch1/sd/elp25/cosmosis:/opt/cosmosis --volume /global/cscratch1/sd/elp25/FireCrown/:/opt/firecrown --volume /global/cscratch1/sd/elp25/TXPipe:/opt/txpipe firecrown run-emcee /global/cscratch1/sd/elp25/FireCrown/examples/des_y1_3x2pt/des_y1_3x2pt.yaml