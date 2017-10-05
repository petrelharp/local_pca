#!/bin/bash
#SBATCH -p long
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH -t 32:00:00
#SBATCH --mem=32G

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

if [ -z ${TAG+x} ]; 
then 
    TAG=$(printf "%06d" $RANDOM);
else
    echo "using existing tag $TAG"
fi

echo "task number $SLURM_ARRAY_TASK_ID and tag $TAG"
echo "looking for tree files in directory $PARAMS"

BASEDIR=$PARAMS
OUTDIR=$BASEDIR/step_mut_$TAG

date
echo "Writing out to $OUTDIR"

ALL_TREEFILES=($BASEDIR/*trees)
ALL_MUTRATES=($(seq .00000001 $(echo "scale=9; .00000009/(${#ALL_TREEFILES[@]}-1)" | bc) .0000001))
TREEFILES=()
MUTRATES=()
for k in $(seq 0 $((${#ALL_TREEFILES[@]}-1)) )
do
    TF=${ALL_TREEFILES[$k]}
    MR=${ALL_MUTRATES[$k]}
    echo "Checking if we've already done $TF"
    if [ ! -f $OUTDIR/$(basename ${TF%trees}vcf) ]
    then
        echo "We haven't!  Will do with mut rate $MR"
        TREEFILES+=($TF)
        MUTRATES+=($MR)
    fi
done

mkdir -p $OUTDIR

echo "Doing:"
echo "./msp-add-mutation.py --basedir $OUTDIR --tree_file ${TREEFILES[@]} --njobs 4 --mut_rate ${MUTRATES[@]}"

./msp-add-mutation.py --basedir $OUTDIR --tree_file ${TREEFILES[@]} --njobs 4 --mut_rate ${MUTRATES[@]}

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
echo "Done!"

