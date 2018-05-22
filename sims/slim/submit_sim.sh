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

export PARAMS="-d RECOMBTYPE='$RECOMBTYPE' -d SELTYPE='$SELTYPE' -d OUTDIR='$OUTDIR' -s $TAG"
export OUTDIR="$OUTDIR"
sbatch -o $OUTDIR/run_${TAG}.out -e $OUTDIR/run_${TAG}.out ./run-sim.sbatch 

