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
echo "looking for tree files in directory $PARAMS"

BASEDIR=$PARAMS
OUTDIR=$BASEDIR/flat_mut_$TAG

date
echo "Writing out to $OUTDIR"

mkdir -p $OUTDIR
./msp-add-mutation.py --basedir $OUTDIR --tree_file $BASEDIR/*trees --njobs 1

date
echo "Indexing and converting to BCF"

for VCF in $OUTDIR/*.vcf
do 
    ./bcf_and_index.sh $VCF 
done 

for VCF in $OUTDIR/*.vcf
do
    if [ -e "${VCF%.vcf}.bcf" ]
    then
        rm $VCF
    fi
done

date
echo "Done generating mutations."

LODIR=$OUTDIR/lostruct_1e6_bp_k2; mkdir -p $LODIR
echo "Writing lostruct results to $LODIR"

./run_lostruct.R --input_dir $OUTDIR --sample_info $BASEDIR/samples01.tsv -t bp -s 1e6 -k 2 -o $LODIR


date
echo "Done!"

