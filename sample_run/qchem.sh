#!/bin/bash

source /opt/qcsetup_5.4omp.sh

export QCSCRATCH=/home/scr/$USER

if [ -d "$QCSCRATCH" ]; then
    mkdir -rf "$QCSCRATCH"
fi

qchem -nt "$SLURM_CPUS_PER_TASK" -save "$1" "$2" "$QCSCRATCH"/"$1".j%

rm -rf "$QCSCRATCH"/"$1".j%
