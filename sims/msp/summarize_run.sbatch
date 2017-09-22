#!/bin/bash
#SBATCH -p long
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -t 4:00:00
#SBATCH --mem=100M

module load R/3.4.1 bcftools

cd $SLURM_SUBMIT_DIR

: ${PARAMS?Must define PARAMS}
if [ -z "$PARAMS" ]
then
    echo "Must define PARAMS (is empty)."
    exit 1
fi

TAG=$(printf "%06d" $RANDOM);

echo "task number $SLURM_ARRAY_TASK_ID and tag $TAG"
echo "looking for BCF files in directory $PARAMS"

OUTDIR=$PARAMS

date
echo "Writing lostruct results to $OUTDIR"
echo "Running:"

Rscript -e "templater::render_template('summarize_run.Rmd', change.rootdir=TRUE, output='${OUTDIR}/run_summary.html')"


date
echo "Done!"
