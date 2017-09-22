#!/bin/bash
#SBATCH -p long
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH -t 24:00:00
#SBATCH --mem=8G

module load gcc
module load python3
module load bright/7.3  # for hdf5_18
module load hdf5_18
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

BASEDIR=$PARAMS

WINTYPE=bp
WINLEN=1e6
KPC=2
MISSING=0.0
LODIR=$BASEDIR/lostruct_${WINLEN}_${WINTYPE}_k${KPC}_m_${MISSING}; mkdir -p $LODIR
echo "Writing lostruct results to $LODIR"
echo "Running:"
echo "  ./run_lostruct.R --input_dir $BASEDIR --sample_info $BASEDIR/../samples01.tsv -t $WINTYPE -s $WINLEN -k $KPC -o $LODIR --missing $MISSING --subsample 0.1"

./run_lostruct.R --input_dir $BASEDIR --sample_info $BASEDIR/../samples01.tsv -t $WINTYPE -s $WINLEN -k $KPC -o $LODIR --missing $MISSING --subsample 0.1

date
echo "Done!"

