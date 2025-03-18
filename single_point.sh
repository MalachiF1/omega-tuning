#!/bin/bash

#SBATCH --job-name=omega_tuning
#SBATCH --output=scratch/slurm/slurm_%j.out

#SBATCH --cpus-per-task=16
#SBATCH --threads-per-core=2
#SBATCH --ntasks-per-core=2
#SBATCH --mem=100000

#SBATCH -A tamars-account
#SBATCH -p tamarsq

source /opt/qcsetup_5.4omp.sh

export QCSCRATCH=/home/scr/$USER

if [ -d "$QCSCRATCH" ]; then
    mkdir -rf "$QCSCRATCH"
fi

export OUTPUT_DIR="$(pwd)/scratch/output"
export INPUT_DIR="$(pwd)/scratch/input"

OMEGA="$1"

if [ ! -d "$OUTPUT_DIR" ]; then
    mkdir "$OUTPUT_DIR"
fi

echo "running job $INPUT_DIR/N_$OMEGA.in to $OUTPUT_DIR/N_$OMEGA.out"
echo "running job $INPUT_DIR/P_$OMEGA.in to $OUTPUT_DIR/P_$OMEGA.out"

qchem -nt "$SLURM_CPUS_PER_TASK" -save "$INPUT_DIR"/N_"$OMEGA".in "$OUTPUT_DIR"/N_"$OMEGA".out "$QCSCRATCH"/N_"$OMEGA"_ID"$SLURM_JOG_ID".j%
qchem -nt "$SLURM_CPUS_PER_TASK" -save "$INPUT_DIR"/P_"$OMEGA".in "$OUTPUT_DIR"/P_"$OMEGA".out "$QCSCRATCH"/P_"$OMEGA"_ID"$SLURM_JOG_ID".j%

# cp -r $QCSCRATCH/$filename.j% $OUTPUT_DIR
rm -rf "$QCSCRATCH"/N_"$OMEGA"_ID"$SLURM_JOG_ID".j%
rm -rf "$QCSCRATCH"/P_"$OMEGA"_ID"$SLURM_JOG_ID".j%
