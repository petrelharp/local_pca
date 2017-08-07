#!/bin/bash

NAME="local_modest"
TAG=$(printf "%06d" $RANDOM); 
DIR="${NAME}_${TAG}"
mkdir -p $DIR
echo "Directory: $DIR"
BASE_PARAMS="-n 1 -j 1 -k 2000 -w 10 -m 1e-3 -u 0.0 -o $DIR -N 1000 -m 1e-3"
NPARAMS=12
echo "$NPARAMS chromosomes."
for CHROM_NUM in $(seq $NPARAMS)
do
    echo "Chrom $CHROM_NUM"
    export PARAMS="$BASE_PARAMS -C $CHROM_NUM"
    sbatch -o $DIR/run_${CHROM_NUM}.out -e $DIR/run_${CHROM_NUM}.out ./run-sim.sbatch 
done
