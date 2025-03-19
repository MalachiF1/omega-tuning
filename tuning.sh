#!/bin/bash

#SBATCH --job-name=omega_tuning
#SBATCH --output=slurm_%j.out

#SBATCH --cpus-per-task=16
#SBATCH --threads-per-core=2
#SBATCH --ntasks-per-core=2
#SBATCH --mem=100000

#SBATCH -A tamars-account
#SBATCH -p tamarsq

# source /opt/qcsetup_5.4omp.sh

# this script runs the tuning python script on the cluster

python tuning.py "$@"
