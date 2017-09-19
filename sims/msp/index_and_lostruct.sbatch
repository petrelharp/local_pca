#!/bin/bash
#SBATCH -p long
#SBATCH -n 1
#SBATCH -c 12
#SBATCH -t 6:00:00
#SBATCH --mem-per-cpu=400M
#SBATCH --ntasks-per-core=1


module load gcc
module load python3
module load bright/7.3  # for hdf5_18
module load hdf5_18
module load python3 R/3.4.1 bcftools

cd $SLURM_SUBMIT_DIR

: ${PARAMS?Must define PARAMS}
if [ -z "$PARAMS" ]
then
    echo "Must define PARAMS (is empty)."
    exit 1
fi

echo "task number $SLURM_ARRAY_TASK_ID"
echo "in directory $PARAMS"

BASEDIR=PARAMS
OUTDIR=$BASEDIR/flat_mut
mkdir -p $OUTDIR
./msp-add-mutation.py --basedir $OUTDIR --tree_file $BASEDIR/*trees --njobs 12
(for VCF in $OUTDIR; do ./bcf_and_index.sh $VCF & done); wait
LODIR=$OUTDIR/lostruct_1e6_bp_k2; mkdir -p $LODIR
./run_lostruct.R --input_dir $OUTDIR --sample_info $BASEDIR/samples01.tsv -t bp -s 1e6 -k 2 -o $LODIR


echo "Done!"

