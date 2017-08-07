#!/bin/bash

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

echo "$PARAMS"

/usr/bin/time --format='elapsed: %E / kernel: %S / user: %U / mem: %M' \
     ./msp-sim.py $PARAMS | grep -v REJECTING

echo "Done!"

