#!/bin/bash
#SBATCH -p long
#SBATCH -n 1
#SBATCH -t 72:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --ntasks-per-core=1

module use /projects/apps/shared/modulefiles
module load python3 tskit SLiM

cd $SLURM_SUBMIT_DIR

: ${OUTDIR?Must define OUTDIR}
if [ -z "$OUTDIR" ]
then
    echo "Must define OUTDIR (is empty)."
    exit 1
fi

echo "Running:"
echo "   slim sims.slim"

/usr/bin/time --format='elapsed: %E / kernel: %S / user: %U / mem: %M' \
    slim -d "OUTDIR='$OUTDIR'" sims.slim &> $OUTDIR/slim_run.log


echo "Done!"

