#!/usr/bin/env Rscript
library(optparse)

invocation <- commandArgs()

usage <- "\
\
Does a local PCA analysis.
To run an analysis on windows of 1000 SNPs, for instance: \
    Rscript run_lostruct.R -t snp -s 1000 \
while to run it on windows of 10000 bps: \
    Rscript run_lostruct.R -t bp -s 10000 \
\
Do
    Rscript run_lostruct.R -h\
for options.\
\
"

option_list <- list(
    # input/output
        make_option( c("-i","--input_dir"),   type="character",        help="Directory with input: indexed .bcf or .vcf.gz files. (REQUIRED)"),
        make_option( c("-I","--sample_info"),   type="character",        help="File with columns labeled 'ID' and 'population' describing populations of samples in bcf file.)"),
        make_option( c("-t","--type"),   type="character",             help="Window by SNP or by bp? (REQUIRED)"),
        make_option( c("-s","--size"),   type="integer",               help="Size of the window, in units of type. (REQUIRED)"),
        make_option( c("-k","--npc"),   type="integer",   default=2L,  help="Number of principal components to compute for each window. [default: %default]"),
        make_option( c("-w","--weightfile"), type="character",         help="File name containing weights to apply to PCA."),
        make_option( c("-m","--nmds"),   type="integer",   default=2L, help="Number of principal coordinates (MDS variables) to compute. [default: %default]"),
        make_option( c("-M","--missing"),   type="double", default=0.0, help="Percent data to introduce as missing.  [default: %default]"),
        make_option( c("-S","--subsample"), type="double", default=1.0, help="Percent of individuals to retain, uniformly across populations.  [default: %default]"),
        make_option( c("-o","--outdir"), type="character",             help="Directory to save results to.  [default: lostruct_results/type_%type_size_%size_jobid_%jobid/]"),
        make_option( c("-j","--jobid"),  type="character", default=formatC(1e6*runif(1),width=6,format="d",flag="0"),   help="Unique job id. [default random]")
    )
opt <- parse_args(OptionParser(option_list=option_list,description=usage))
if (is.null(opt$outdir)) { opt$outdir <- file.path("lostruct_results", 
                                   sprintf( "type_%s_size_%d_weights_%s_jobid_%s", 
                                           opt$type, 
                                           opt$size, 
                                           if (is.null(opt$weightfile)) { "none" } else { gsub("[.].*","",basename(opt$weightfile)) }, 
                                           opt$jobid ) ) }
if (is.null(opt$input_dir) || is.null(opt$type) || is.null(opt$size)) { stop(usage) }

opt$start.time <- Sys.time()
opt$run.dir <- normalizePath(".")

if (is.null(opt$sample_info)) {
    warning("No sample information file - report from summarize_run.Rmd won't have colors.")
} else if  (!file.exists(opt$sample_info)) {
    stop(sprintf("Sample info file %s does not exist.",opt$sample_info))
} else {
    opt$sample_info <- normalizePath(opt$sample_info)
}

# subsamples
if (opt$subsample < 1.0) {
    if (is.null(opt$sample_info)) { stop("Cannot subsample without a sample file.") }
    samps <- read.table(opt$sample_info, sep="\t", header=TRUE)
    names(samps) <- tolower(names(samps))
    # hack for msprime output
    if (is.numeric(samps$id)) { samps$id <- factor(paste0("msp_", samps$id)) }
    subsamps <- as.character(samps$id)[sort(unlist(tapply(1:nrow(samps), samps$population, 
                    function (x) sample(x, ceiling(opt$subsample*length(x))))))]
}

# weights
if (is.null(opt$weightfile)) {
    opt$weights <- 1
} else {
    if (!file.exists(opt$weightfile)){
        stop(sprintf("Weight file %s does not exist.",opt$weightfile))
    }
    opt$weights <- read.table(opt$weightfile,sep='\t',header=TRUE)$weight
}

# VCF files
if (!dir.exists(opt$input_dir)) {
    stop(sprintf("Input directory %s does not exist.",opt$input_dir))
}
bcf.files <- list.files(opt$input_dir,".*.(bcf|vcf.gz|vcf.bgz)$",full.names=TRUE)
if (length(bcf.files)==0) {
    stop(sprintf("No bcf or vcf.gz files found in input directory %s",opt$input_dir))
}
bcf.files <- normalizePath(bcf.files)
names(bcf.files) <- make.names(gsub("[.](bcf|vcf.gz|vcf.bgz)$","",basename(bcf.files)))

dir.create( opt$outdir, showWarnings=FALSE, recursive=TRUE )

# override vcf_windower to introduce missing data if desired
if (is.numeric(opt$missing) && (opt$missing > 0)) {
    vcf_windower <- function (...) {
        f <- lostruct::vcf_windower(...)
        g <- as.winfun(f=function (...) {
                            out <- f(...);
                            m <- (rbinom(length(out), size=1, prob=opt$missing) > 0);
                            out[m] <- NA; return(out) },
                       max.n=attr(f, "max.n"),
                       samples=attr(f, "samples"),
                       region=attr(f, "region"))
        return(g)
    }
}

# setup
library(lostruct)
options(datatable.fread.input.cmd.message=FALSE)

all.pcas <- numeric(0)       # will be a numeric matrix of eigen values/vectors
all.lengths <- numeric(0)    # will be a numeric vector of numbers of windows per chromosome
all.regions <- data.frame()  # will be a data.frame of the chromsome, start, stop for each window
chroms <- c()                # a list of chromosome names
all.files <- c()             # a list of the files the chromosomes appear in

# local PCA, by chromosome
for (k in seq_along(bcf.files)) {
    bcf.file <- bcf.files[k]

    file.positions <- vcf_positions(bcf.file)
    for (chrom.index in seq_along(file.positions)) {
        chrom.name <- names(file.positions)[chrom.index]
        # prefix with 'chr' if numeric
        if (suppressWarnings(!is.na(as.numeric(chrom.name)))) {
            chrom.name <- paste0('chr', chrom.name)
        }
        if (chrom.name %in% chroms) {
            stop(sprintf("Chromosome %s appears in more than one file.", chrom.name))
        }
        chroms <- c(chroms, chrom.name)
        all.files <- c(all.files, bcf.file)
        pca.file <- file.path( opt$outdir, paste0(chrom.name,".pca.csv") )
        regions.file <- file.path( opt$outdir, paste0(chrom.name,".regions.csv") )
        if (file.exists(pca.file)) { 
            warning(paste("File",pca.file,"already exists! Not recomputing.")) 
            pca.stuff <- as.matrix( data.table::fread(pca.file,header=TRUE) )
            these.regions <- data.table::fread(regions.file,header=TRUE)
        } else {
            cat("Finding PCs for chromsome", chrom.name,
                "in file", bcf.file, "and writing out to", pca.file, "and", regions.file, "\n")
            sites <- file.positions[chrom.index]
            if (opt$subsample < 1.0) {
                win.fn <- vcf_windower(bcf.file, size=opt$size, type=tolower(opt$type),
                                       samples=subsamps, sites=sites)
            } else {
                win.fn <- vcf_windower(bcf.file, size=opt$size, type=tolower(opt$type),
                                       sites=sites)
            } 
            these.regions <- region(win.fn)()
            system.time( 
                        pca.stuff <- eigen_windows( win.fn, k=opt$npc, w=opt$weights ) 
                    )
            write.csv( pca.stuff, file=pca.file, row.names=FALSE )
            write.csv( these.regions, file=regions.file, row.names=FALSE )
        }
        all.pcas <- rbind( all.pcas, pca.stuff )
        all.lengths <- c(all.lengths, nrow(pca.stuff))
        all.regions <- rbind( all.regions, these.regions )
    }
}
names(all.lengths) <- chroms
rm(pca.stuff)

# write out the config json
opt$bcf_files <- all.files
opt$chrom.names <- chroms
cat( jsonlite::toJSON( opt, pretty=TRUE ), file=file.path( opt$outdir, "config.json" ) )

# distance matrix
cat("Done finding PCs, computing distances.\n")
system.time( pc.distmat <- pc_dist( all.pcas, npc=opt$npc ) )

# MDS on the resulting distance matrix
mds.file <- file.path( opt$outdir, "mds_coords.csv" )
cat("Done computing distances, running MDS and writing results to", mds.file, "\n")
na.inds <- is.na( all.pcas[,1] ) # there may be windows with missing data
mds.coords <- cbind( data.frame( 
                        chrom=rep(chroms,all.lengths),
                        window=unlist(lapply(all.lengths,seq_len)) ),
                        cmdscale( pc.distmat[!na.inds,!na.inds], k=opt$nmds )[ ifelse( na.inds, NA, cumsum(!na.inds) ), ]
                    )
colnames(mds.coords)[-(1:2)] <- paste0("MDS",seq_len(opt$nmds))
write.csv( mds.coords, mds.file, row.names=FALSE )

cat("All done!")
Sys.time()
