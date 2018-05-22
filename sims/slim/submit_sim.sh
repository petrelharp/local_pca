#!/bin/bash

USAGE="
Usage:
   $0 SELTYPE RECOMBTYPE
"

if [ $# -lt 2 ]
then
    echo "$USAGE"
    exit 0
fi

SELTYPE="$1"
RECOMBTYPE="$2"

NAME="run"
TAG=$(printf "%06d" $RANDOM); 
OUTDIR="${NAME}_${TAG}"
mkdir -p $OUTDIR
echo "Directory: $OUTDIR"

echo "defineConstant(\"RECOMBTYPE\", \"${RECOMBTYPE}\");" >> $OUTDIR/parameters.slim
echo "defineConstant(\"SELTYPE\", \"${SELTYPE}\");" >> $OUTDIR/parameters.slim

export OUTDIR
sbatch -o $OUTDIR/slurm_${TAG}.out -e $OUTDIR/slurm_${TAG}.out ./run_sim.sbatch 

