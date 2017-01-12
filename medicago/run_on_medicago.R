#!/usr/bin/Rscript --vanilla
library(optparse)

invocation <- commandArgs()

usage <- "\
\
Does a local PCA analysis on the Medicago truncatula hapmap data (see medicago_data_setup.html for details). \
To run an analysis on windows of 1000 SNPs, for instance: \
    Rscript run_on_medicago.R -t snp -s 1000 \
while to run it on windows of 1000 bps: \
    Rscript run_on_medicago.R -t bp -s 1000 \
\
Do
    Rscript run_on_medicago.R -h\
for options.\
\
"

option_list <- list(
    # input/output
        make_option( c("-t","--type"),   type="character",             help="Window by SNP or by bp?"),
        make_option( c("-s","--size"),   type="integer",               help="Size of the window, in units of type."),
        make_option( c("-k","--npc"),   type="integer",   default=2L,  help="Number of principal components to compute for each window. [default: %default]"),
        make_option( c("-w","--weightfile"), type="character",            help="File name containing weights to apply to PCA."),
        make_option( c("-m","--nmds"),   type="integer",   default=2L, help="Number of principal coordinates (MDS variables) to compute. [default: %default]"),
        make_option( c("-o","--outdir"), type="character",             help="Directory to save results to.  [default: lostruct/results_type_%type_size_%size_jobid_%jobid/]"),
        make_option( c("-j","--jobid"),  type="character", default=formatC(1e6*runif(1),width=6,format="d",flag="0"),   help="Unique job id. [default random]")
    )
opt <- parse_args(OptionParser(option_list=option_list,description=usage))
if (is.null(opt$outdir)) { opt$outdir <- file.path("lostruct", 
                                   sprintf( "results_type_%s_size_%d_weights_%s_jobid_%s", 
                                           opt$type, 
                                           opt$size, 
                                           if (is.null(opt$weightfile)) { "none" } else { gsub("[.].*","",basename(opt$weightfile)) }, 
                                           opt$jobid ) ) }
if (is.null(opt$type) || is.null(opt$size)) { stop(usage) }

opt$start.time <- Sys.time()
opt$run.dir <- normalizePath(".")
print(opt) # this will go in the logfile

# weights
if (is.null(opt$weightfile)) {
    opt$weights <- 1
} else {
    opt$weights <- read.table(opt$weightfile,sep='\t',header=TRUE)$weight
}

dir.create( opt$outdir, showWarnings=FALSE, recursive=TRUE )
cat( jsonlite::toJSON( opt, pretty=TRUE ), file=file.path( opt$outdir, "config.json" ) )

# setup
library(lostruct)

# VCF files
chroms <- paste0("chr",1:8)
bcf.files <- file.path( "data", paste0(chroms,"-filtered-set-2014Apr15.bcf") )
names(bcf.files) <- chroms

all.pcas <- numeric(0)       # will be a numeric matrix of eigen values/vectors
all.lengths <- numeric(0)    # will be a numeric vector of numbers of windows per chromosome
all.regions <- data.frame()  # will be a data.frame of the chromsome, start, stop for each window

# local PCA, by chromosome
for (bcf.file in bcf.files) {
    pca.file <- file.path( opt$outdir, sprintf(gsub(".bcf",".pca.csv",basename(bcf.file))) )
    regions.file <- file.path( opt$outdir, sprintf(gsub(".bcf",".regions.csv",basename(bcf.file))) )
    if (file.exists(pca.file)) { 
        warning(paste("File",pca.file,"already exists! Not recomputing.")) 
        pca.stuff <- as.matrix( data.table::fread(pca.file,header=TRUE) )
        these.regions <- data.table::fread(regions.file,header=TRUE)
    } else {
        cat("Finding PCs for", bcf.file, "and writing out to", pca.file, "and", regions.file, "\n")
        win.fn <- vcf_windower(bcf.file, size=opt$size, type=tolower(opt$type) )
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
names(all.lengths) <- chroms
rm(pca.stuff)

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
