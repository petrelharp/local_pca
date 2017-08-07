#!/bin/bash
#SBATCH -n 12
#SBATCH -t 48:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --ntasks-per-core=1
#SBATCH -J background_modest

NAME="background_modest"
TAG=$(printf "%06d" $RANDOM); 
DIR="${NAME}_${TAG}"
BASE_PARAMS="-n 1 -j 1 -k 2000 -w 10 -m 1e-3 -u 0.0 -o $DIR"
CHROM_NUM=0
for N in 200 200 400 400 600 600 800 800 1000 1000 1200 1200
do
    SEED=$(printf "%06d" $RANDOM); 
    echo "Chrom number $CHROM_NUM - seed $SEED"
    PARAMS="$BASE_PARAMS -d $SEED -N $N -C $CHROM_NUM"
    srun ./run-sim.sh &>  $DIR/run.log
    ((CHROM_NUM++))
done
