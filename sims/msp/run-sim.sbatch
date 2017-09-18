#!/bin/bash
#SBATCH -p long
#SBATCH -n 1
#SBATCH -t 72:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --ntasks-per-core=1


module load gcc
module load python3
module load bright/7.3  # for hdf5_18
module load hdf5_18

cd $SLURM_SUBMIT_DIR

: ${PARAMS?Must define PARAMS}
if [ -z "$PARAMS" ]
then
    echo "Must define PARAMS (is empty)."
    exit 1
fi

SEED=$(printf "%06d" $RANDOM); 

echo "Chrom number $SLURM_ARRAY_TASK_ID - seed $SEED"
ALL_PARAMS="$PARAMS -d $SEED"
echo "$ALL_PARAMS"

/usr/bin/time --format='elapsed: %E / kernel: %S / user: %U / mem: %M' \
     ./msp-sim.py $ALL_PARAMS | grep -v REJECTING

echo "Done!"

