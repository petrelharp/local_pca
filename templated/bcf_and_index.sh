#!/bin/bash

if [ $# -lt 1 ]
then
    echo "Usage:"
    echo "   $0 [vcf file [vcf file [...]]]"
    echo "For each file on the command line, will convert to bcf and index it."
fi

while (( "$#" )) 
do
    if [[ ! -e "$1" ]] 
    then 
        echo "File $1 does not exist."
    else 
        OUT=${1%.vcf}.bcf
        echo "Converting $1 to $OUT (and indexing)."
        bcftools convert -o $OUT -O b $1
        bcftools index $OUT
    fi
    shift
done
