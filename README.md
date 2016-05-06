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

There are standalone examples for each of the three datasets studied:

### [POPRES](popres/) (*Homo sapiens*, SNP chip data from a few worldwide populations)

Chromosome 1 is the example given .

### [DPGP](dpgp/) (*Drosophila melanogaster* population genome project)

Chromosome 3L is the example given .

### [Medicago](medicago/) (*Medicago truncatula* hapmap)

For Medicago, it calculates the pairwise distance for all 8 chromosome together and then apply MDS and use subset of the whole MDS result for each chromosome. 


# Using the R package

The example scripts in the directories above mostly work *without* the R package.
To start using the code on your own data, have a look at these files:

* [Setting up the medicago data](medicago/medicago_data_setup.html) : after documenting where the data are from,
    does local PCA on a small subset of the whole dataset, to establish how the functions work.

* [Script for medicago analysis](medicago/run_on_medicago.R) : an Rscript to run the same analysis on medicago data,
    varying various parameters by command-line options: run `Rscript run_on_medicago.R --help` for a list.
