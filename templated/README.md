# Templated reports with lostruct

There are a few parameters to play around with in lostruct.
Here's an example of how to do this in an easy way,
and compare between parameter combinations.

In this directory:

* `run_lostruct.R` - an R script that computes everything we need.
* `summarize_run.Rmd` - an Rmarkdown file that can be compiled to visualize the results.

## Overview

The way this works is:

1. Put your data, as indexed bcf files, into the `data/` directory.

2. Run `run_lostruct.R` to

    a. Find local PCs in windows, separately for each bcf file.
    b. Compute the distance matrix between windows.
    c. Run MDS on that distance matrix

    The results will be output to a directory with a unique name, 
    that includes a `.json` file with the configuration parameters in it.

    `./run_lostruct.R -h` will give you a list of command-line options.

3. Put information about your samples into `data/sample_info.tsv` - tab-separated, with columns `ID` and `population`.

4. Compile `summarize_run.Rmd` in that directory to produce an html report with pretty figures.

## Prerequisites

To compile the report you'll probably want [templater](https://github.com/petrelharp/templater).
In R:
```r
library(devtools)
install_github("petrelharp/templater")
```

## Example

Using windows of 20 adjacent SNPs
(note: these windows are far too small for real data),
compute everything:

```bash
./run_lostruct.R -i data -t snp -s 20 -I data/sample_info.tsv -j 0001
# ...
# Finding PCs for data/chr1.bcf and writing out to lostruct_results/type_snp_size_20_weights_none_jobid_0001/chr1.pca.csv and lostruct_results/type_snp_size_20_weights_none_jobid_0001/chr1.regions.csv 
# Finding PCs for data/chr2.bcf and writing out to lostruct_results/type_snp_size_20_weights_none_jobid_0001/chr2.pca.csv and lostruct_results/type_snp_size_20_weights_none_jobid_0001/chr2.regions.csv 
# Warning messages:
# 1: In vcf_windower_snp(file = file, sites = sites, size = size, samples = samples) :
#   Trimming from chromosome ends: 1: 14 SNPs.
# 2: In vcf_windower_snp(file = file, sites = sites, size = size, samples = samples) :
#   Trimming from chromosome ends: 1: 14 SNPs.
# Done finding PCs, computing distances.
#    user  system elapsed 
#   0.268   0.012   0.280 
# Done computing distances, running MDS and writing results to lostruct_results/type_snp_size_20_weights_none_jobid_0001/mds_coords.csv 

```

*(Note: we're setting the jobid explicitly here, with `-j 0001`, but recommend leaving it unset, so that it is randomly generated, and you don't overwrite previous results.)*
This may take a long time on your real dataset.
The warning is because since the number of SNPs isn't evenly divisible by 20, it discards the dangling 14 at the end of each chromosome.
It's put everything in the directory `lostruct_results/type_snp_size_20_weights_none_jobid_0001/`:

```bash
ls lostruct_results/type_snp_size_20_weights_none_jobid_0001/
# chr1.pca.csv  chr1.regions.csv  chr2.pca.csv  chr2.regions.csv  config.json  mds_coords.csv
```
These files are:

* `config.json` - the options passed to `run_lostruct.R`, so you know what you did (also used in making the report).
* `chr1.pca.csv` - the percent variation explained and top two eigenvalues and eigenvectors of the windows, one window per row - output by `eigen_windows()`.
* `chr1.regions.csv` - the locations of the windows, one per row: chromosome, first SNP position, last SNP position.
* `mds_coords.csv` - the MDS coordinates of the windows - file name, window number, and MDS coordinates.

The `589131` suffix on the directory is to make the directory unique - override this with the `-o` option.

Now we want to compile the template.  You can do this from the command line as
```bash
Rscript -e 'templater::render_template("summarize_run.Rmd",output="lostruct_results/type_snp_size_20_weights_none_jobid_0001/run_summary.html",change.rootdir=TRUE)'
```
This will produce the file [lostruct_results/type_snp_size_20_weights_none_jobid_0001/run_summary.html](lostruct_results/type_snp_size_20_weights_none_jobid_0001/run_summary.html),
which you can view in a web browser.

If you want to examing things in more detail,
open up R (in this directory, not in the results directory!)
and run
```r
library(templater)
render_template(
        "summarize_run.Rmd",
        output="lostruct_results/type_snp_size_20_weights_none_jobid_0001/run_summary.html",
        change.rootdir=TRUE,
        envir=environment())
```
Now, you'll be able to see everything that is computed while parsing the Rmd file.
This makes it easy for you to add things to the report
or modify it to output figures as pdf for publication.

## Want pdfs?

The report will make pdf versions of the figures, also -- just edit `summarize_run.Rmd`
and change `do.pdfs <- FALSE` to `do.pdfs <- TRUE`.
The pdfs will show up in the figure subdirectory for the report
(also listed in the report under the figure if you turn this on.
In the example above it is `lostruct_results/type_snp_size_20_weights_none_jobid_0001/figure/run_summary/`.
