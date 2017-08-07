#!/bin/bash

NAME="background_strong"
TAG=$(printf "%06d" $RANDOM); 
DIR="${NAME}_${TAG}"
mkdir -p $DIR
echo "Directory: $DIR"
BASE_PARAMS="-n 1 -j 1 -k 2000 -w 10 -m 1e-3 -u 0.0 -o $DIR"
MORE_PARAMS=("-N 200" "-N 200" "-N 400" "-N 400" "-N 600" "-N 600" "-N 800" "-N 800" "-N 1000" "-N 1000" "-N 1200" "-N 1200")
NPARAMS=${#MORE_PARAMS[@]}
echo "$NPARAMS chromosomes."
for CHROM_NUM in $(seq $((NPARAMS)))
do
    echo "Chrom $CHROM_NUM"
    export PARAMS="$BASE_PARAMS ${MORE_PARAMS[$((CHROM_NUM-1))]} -C $CHROM_NUM"
    sbatch -o $DIR/run_${CHROM_NUM}.out -e $DIR/run_${CHROM_NUM}.out ./run-sim.sbatch 
done
