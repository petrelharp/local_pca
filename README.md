Local PCA/population structure (lostruct)
=========================================

To install the package, make sure you have `devtools` (by doing `install.packages("devtools")`),
and then running 
```
install.packages("data.table")
devtools::install_github("petrelharp/local_pca/lostruct")
library(lostruct)
```
*Note:* the library is called `lostruct`.

## Using the R package

The example scripts in the directories above mostly work *without* the R package.
To start using the code on your own data, have a look at these files:

* [A quick example](popres/popres_example.R) : in four lines of code, reads in chromosome 22 from a TPED, and does local PCA.

* [Setting up the medicago data](medicago/medicago_data_setup.html) : after documenting where the data are from,
    does local PCA on a small subset of the whole dataset, to establish how the functions work.

* [Script for medicago analysis](medicago/run_on_medicago.R) : an Rscript to run the same analysis on medicago data,
    varying various parameters by command-line options: run `Rscript run_on_medicago.R --help` for a list.

* [Report summarizing an analysis](medicago/summarize_run.Rmd) : an Rmarkdown file that can be compiled with [templater](https://github.com/petrelharp/templater) to produce visualizations of the results of the above.

## Prerequisites

- To use the functions to read in windows out of BCF file,
    you will need [bcftools](http://www.htslib.org/doc/bcftools.html).
- To compile the example report, you probably want 
    [templater](https://github.com/petrelharp/templater).

## Standalone code

Also included is code we used to analyze the datasets in the paper (before the R package was written).
The general order to see the code in each directory is 

1. recode : turn bases into numbers
2. PCA : find local PCs
3. distance : compute distance matrix between windows from local PCs
4. MDS : visualize the result

There are standalone examples for each of the three datasets studied:

### [POPRES](popres/) (*Homo sapiens*, SNP chip data from a few worldwide populations)

Chromosome 1 is the example given.  See also [popres_example.R](popres/popres_example.R) for an example of some steps using the package.

- [POPRES_SNPdata_recode12.R](popres/POPRES_SNPdata_recode12.R) : recodes the TPED as numeric
- [POPRES_cov.R](popres/POPRES_cov.R) : computes covariance matrix for the entire chromsome 1
- [POPRES_PCA_win100.R](popres/POPRES_PCA_win100.R) : computes local PCs
- [POPRES_jackknife_var.R](popres/POPRES_jackknife_var.R) : estimates SE of local PCs
- [POPRES_distance.R](popres/POPRES_distance.R) : computes distance matrix from local PCs
- [POPRES_MDS.R](popres/POPRES_MDS.R) : finds and plots MDS visualization of distance matrix

### [DPGP](dpgp/) (*Drosophila melanogaster* population genome project)

Chromosome 3L is the example given .

- [DPGP_recode_and_cov.R](dpgp/DPGP_recode_and_cov.R) : recodes data as numeric, removes individuals with more than 8% missing data, sites with more than 20% missing data, and computes whole-chromosome covariance matrix
- [DPGP_PCA_plot.R](dpgp/DPGP_PCA_plot.R) : plots PCs for entire 3L
- [DPGP_PCA_win103.R](dpgp/DPGP_PCA_win103.R) : computes local PCs along 3L in windows of 1000 SNPs
- [DPGP_var_between_win.R](dpgp/DPGP_var_between_win.R) : computes variance of PCs between windows
- [DPGP_jackknife_var.R](dpgp/DPGP_jackknife_var.R) : does jackknife estimate of SE for PCs on windows of 1000 SNPs
- [DPGP_distance.R](dpgp/DPGP_distance.R) : computes distance matrix from local PCs
- [DPGP_MDS_1d.R](dpgp/DPGP_MDS_1d.R) : computes and plots MDS plots from distance matrix
- [DPGP_get_extreme_points.R](dpgp/DPGP_get_extreme_points.R) : identifies extreme points (with interaction)
- [DPGP_combine_extremes_and_get_cov.R](dpgp/DPGP_combine_extremes_and_get_cov.R) : combines each of three sets of extreme windows and computes covariances for each

### [Medicago](medicago/) (*Medicago truncatula* hapmap)

For Medicago, it calculates the pairwise distance for all 8 chromosome together and then apply MDS and use subset of the whole MDS result for each chromosome. 

- [Medicago_VCF_recode.py](medicago/Medicago_VCF_recode.py) : recodes VCF file as numeric
- [Medicago_recode_and_cov.R](medicago/Medicago_recode_and_cov.R) : computes covariance matrix for (entire) chromosome 1
- [Medicago_PCA_win104.R](medicago/Medicago_PCA_win104.R) : computes local PCs for chromosome 1
- [Medicago_distance_all_chr.R](medicago/Medicago_distance_all_chr.R) : computes a distance matrix from PC information
- [Medicago_MDS.R](medicago/Medicago_MDS.R) : computes and plots MDS plots from the distance matrix


# A note on implementation:

This method works through the genome doing something (PCA on the covariance matrix)
one window at a time.  Because of this, it can be frustratingly slow to first load
the entire dataset into memory.  There are several methods implemented here to avoid this;
for instance, `vcf_windower()` which is used to [compute PCs for the medicago data](medicago/run_on_medicago.R).
The interface is via a function that takes an integer, `n`,
and returns a data frame of the genomic data in the `n`th window.

