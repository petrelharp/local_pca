Local PCA/population structure (lostruct)
=========================================

To load the package, make sure you have `devtools` (by doing `install.packages("devtools")`),
and then running `library(devtools); load_all("PATH/TO/THIS/DIRECTORY/package")`.

## Workflow

The general order to see the code is 

1. recode
2. cov
3. PCA
4. distance
5. MDS

## POPRES (human, SNP chip)

Chromosome 1 is the example given .

## DPGP (Drosophila)

Chromosome 3L is the example given .

## Medicago

For Medicago, it calculates the pairwise distance for all 8 chromosome together and then apply MDS and use subset of the whole MDS result for each chromosome. 
