#!/bin/bash

#SBATCH --account kireevlab
#SBATCH --job-name standardize
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem-per-cpu=10G
#SBATCH --output info.log
#SBATCH --time 4:00:00

python3 standardize.py molport_10k.smi standardize.smi
