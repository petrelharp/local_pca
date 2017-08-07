#!/bin/bash
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -t 24:00:00
#SBATCH --mem-per-cpu=8G

module load gcc
module load python3
module load bright/7.3  # for hdf5_18
module load hdf5_18

cd $SLURM_SUBMIT_DIR

for N in 1000;
do
    for L in 2e6 4e6 5e6 10e6
    do
        X="test_mN4_${L}_${N}_"
        SEED=$(printf "%06d" $RANDOM); 
        mkdir ${X}$SEED; 
        (/usr/bin/time --format='elapsed: %E / kernel: %S / user: %U / mem: %M' \
         ./msp-sim.py -n 1 -k 2000 -w 10 -N $N -o ${X}$SEED -d $SEED -u 0.0 -j 1\
            -m .001 -L $L \
            | grep -v REJECTING ) &> ${X}$SEED/run.log  &
    done
    wait
done

echo "Done!"
