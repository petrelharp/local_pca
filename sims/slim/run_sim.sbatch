#!/bin/bash
#SBATCH -p long
#SBATCH -n 1
#SBATCH -t 72:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --ntasks-per-core=1

module use /projects/apps/shared/modulefiles
module load python3 tskit SLiM

cd $SLURM_SUBMIT_DIR

: ${PARAMS?Must define PARAMS}
if [ -z "$PARAMS" ]
then
    echo "Must define PARAMS (is empty)."
    exit 1
fi

: ${OUTDIR?Must define OUTDIR}
if [ -z "$OUTDIR" ]
then
    echo "Must define OUTDIR (is empty)."
    exit 1
fi

echo "Running:"
echo "   slim $PARAMS sims.slim"

/usr/bin/time --format='elapsed: %E / kernel: %S / user: %U / mem: %M' \
    slim "$PARAMS" sims.slim &> $OUTDIR/run.log


echo "Done!"

