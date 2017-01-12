
#' Make Regions String from Data Frame of Regions
#'
#' Given a data frame with variables \code{chrom}, \code{start}, and \code{end},
#' make a regions string that can be passed to bcftools.
#'
#' @param regions A data frame with \code{chrom}, \code{start}, and \code{end}.
#' @return A character string of the form "chr1:start1-end1,chr2:start2,end2".
#' @export
region_string <- function (regions) {
    if (nrow(regions)==0) { return( NULL ) }
    paste( paste0( regions$chrom, ":", as.numeric(regions$start), "-", as.numeric(regions$end) ), collapse="," )
}


#' Query a VCF File For Numeric Genotypes
#'
#' Uses bcftools to query vcf or bcf files for genotypes in a collection of regions,
#' returning the genotypes in numeric format, giving the number of nonreference alleles.
#' Assumes the data are in the format either "0/1", "0|1" (for diploids), or "0" (for haploids).
#' Data are returned in the order encountered in the VCF file even if the regions are not in this order.
#'
#' @param file The name of the indexed vcf file.
#' @param regions A data frame containing chrom, start, and end of the regions to be extracted.
#' @param samples A character vector of sample names to be extracted.
#' @param verbose Whether to output with \code{cat()} the bcftools calls.
#' @param recode Whether to recode genotypes as numeric (assuming diallelic sites only).
#' @return An integer matrix, giving the number of nonreference alleles for each sample.
#' For instance, "0/0", "0|0", and "0" are coded as 0, while "0/1", "2|0", or "1" are coded as 1,
#' and "1/1", "1/2" or "1|4" are all coded as 2.
#'
#' The format is determiend by looking at the first 100 loci in the first 100 individuals.
#' @export
vcf_query <- function (file, regions, samples, verbose=FALSE, recode=TRUE) {
	bcf.args <- c("bcftools", "query", "-f", "'[ %GT]\\n'")
	if (!missing(regions) && length(regions)>0) { bcf.args <- c( bcf.args, "-r", region_string(regions) ) }
	if (!missing(samples) && length(samples)>0) { bcf.args <- c( bcf.args, "-s", paste(samples,collapse=',') ) }
	bcf.args <- c( bcf.args, file )
    bcf.call <- paste(bcf.args, collapse=" ")
    if (verbose) cat(bcf.call, "\n")
    gt.text <- tryCatch( as.matrix(data.table::fread( bcf.call, header=FALSE, sep=' ', data.table=FALSE )),
                   error=function (e) { if ( grepl("File is empty", e$message) ) { NULL } else { stop(paste("Error. Is bcftools installed?\n",e)) } } )
    if (is.null(gt.text) || !recode) { return(gt.text) }
    if (length(grepl("[0-9]\\|[0-9]", gt.text[seq_len(min(nrow(gt.text),100)),seq_len(min(ncol(gt.text),100))]))>0) {
        gt <- c(0L,1L,1L,2L)[match( unlist(gt.text), c("0|0","0|1","1|0","1|1") )]
        gt[ grep( "([2-9]\\|0)|(0\\|[2-9])", gt.text ) ] <- 1L
        gt[ grep( "([2-9]\\|1)|(1\\|[2-9])|([2-9]\\|[2-9])", gt.text ) ] <- 2L
    } else if (length(grepl("[0-9]\\/[0-9]", gt.text[seq_len(min(nrow(gt.text),100)),seq_len(min(ncol(gt.text),100))]))>0) {
        gt <- c(0L,1L,1L,2L)[match( unlist(gt.text), c("0/0","0/1","1/0","1/1") )]
        gt[ grep( "([2-9]\\/0)|(0\\/[2-9])", gt.text ) ] <- 1L
        gt[ grep( "([2-9]\\/1)|(1\\/[2-9])|([2-9]\\/[2-9])", gt.text ) ] <- 2L
    } else {
        gt <- c(0L,1L)[match( unlist(gt.text), c("0","1") )]
        gt[ grep( "[2-9]", gt.text ) ] <- 1L
    }
    dim(gt) <- dim(gt.text)
    return(gt)
}

#' Query Multiple VCF Files
#'
#' Applies \code{vcf_query} to multiple VCF files: must specify which chromosomes are in which file with \code{chrom.list}.
#'
#' @param chrom.list A list of character vectors of chromosome names, the \code{n}th saying which chromosomes are available in the \code{n}th vcf file.
#' @param file The name(s) of the indexed vcf file(s).
#' @param regions A data frame containing chrom, start, and end of the regions to be extracted.
#' @param ... Other arguments passed to \code{vcf_query}.
#' @return The results of the separate \code{vcf_query} calls, \code{rbind}'ed together
#' (this will *not* necessarily be in the same order as \code{regions}!).
#' @export
multi_vcf_query <- function (chrom.list, file=names(chrom.list), regions, ... ) {
    chroms <- unlist(chrom.list)
    files <- file[ rep(seq_along(chrom.list),sapply(chrom.list,length)) ]
    do.call( rbind, lapply( seq_along(file), function (kf) {
                this.regions <- regions[ files[match(regions$chrom,chroms)]==file[k], ]
                vcf_query( file[k], this.regions, ... )
        } ) )
}

#' Function to Query Multiple VCF Files
#'
#' Returns a function that when called with an integer \code{n} will return the \code{n}th region.
#'
#' @param chrom.list A list of character vectors of chromosome names, the \code{n}th saying which chromosomes are available in the \code{n}th vcf file.
#' @param file The name(s) of the indexed vcf file(s).
#' @param regions A data frame containing chrom, start, and end of the regions to be extracted.
#' @param ... Other arguments passed to \code{vcf_query}.
#' @return A function, say, \code{f}, so that \code{f(n)} is the result of \code{n}th \code{vcf_query} call.
#' @export
multi_vcf_query_fn <- function (chrom.list, file=names(chrom.list), regions, ... ) {
    chroms <- unlist(chrom.list)
    files <- file[ rep(seq_along(chrom.list),sapply(chrom.list,length)) ]
    qfun <- function (n) {
        this.region <- regions[n,]
        vcf_query( files[match(this.region$chrom,chroms)], this.region, ... )
    }
    attr(qfun,"max.n") <- nrow(regions)
    attr(qfun,"region") <- regions
    attr(qfun,"samples") <- samples(files[1])
    return(qfun)
}


#' Read Position From a VCF File
#'
#' Reads the CHROM and POS columns of a VCF file, with bcftools,
#' returning a list of positions named by the corresponding chromosomes.
#'
#' @param file The name of the indexed vcf files.
#' @return A named list of integer vectors; names correspond to chromosomes.
#' @export
vcf_positions <- function (file) {
	# bcf.con <- pipe(paste("bcftools query -f '%CHROM\\t%POS\\n'",file),open="r")
	# bcf.sites <- read.table(bcf.con, sep='\t')
	# close(bcf.con)
	bcf.sites <- data.table::fread(paste("bcftools query -f '%CHROM\\t%POS\\n'",file), header=FALSE, sep='\t', data.table=FALSE)
	colnames(bcf.sites) <- c("chrom","pos")
	return( tapply( bcf.sites$pos, bcf.sites$chrom, identity) )
}

#' Find Sample IDs from a VCF File
#'
#' Returns the character vector of samples (by calling \code{bcftools query -l}).
#'
#' @param file The name of the indexed vcf file.
#' @return A character vector of sample IDs.
#' @export
vcf_samples <- function (file) {
    return( system2( "bcftools", c("query","-l",file), stdout=TRUE ) )
}


#' Construct a Window Extractor for a VCF File
#'
#' Returns a function that, given an integer \code{n}, returns a numeric matrix of genotypes from the the \code{n}th window of the VCF file.
#'
#' @param file The name of the indexed vcf file.
#' @param size The size of the window.
#' @param type The units of the window: 'bp' or 'snp'?
#' @param sites The positions in the VCF file, as returned by \code{vcf_positions}.
#' @param samples A character vector of sample names to be extracted.
#' @return A class "winfun" window extractor function that returns an integer matrix, giving the number of alternate alleles for each sample.
#' Such functions also have three attributes: \code{max.n}, giving the index of the largest window,
#' \code{samples}, the sample IDs that are extracted (corresponding to the columns of the matrix that is returned),
#' and \code{region}, which is a function that takes an integer vector and returns a data frame giving chromosome, start, and end of the corresponding windows.
#' @export
vcf_windower <- function (
			file, 
			size, 
			type, 
			sites=vcf_positions(file),
			samples=vcf_samples(file)) {
	if (type=="bp") {
		vcf_windower_bp( file=file, sites=sites, size=size, samples=samples )
	} else if (type=="snp") {
		vcf_windower_snp( file=file, sites=sites, size=size, samples=samples )
	}
}

vcf_windower_bp <- function (file, sites, size, samples=vcf_samples(file)) {
	chroms <- names(sites)
	chrom.starts <- sapply( sites, min )
	chrom.lens <- sapply( sites, max ) - chrom.starts
	chrom.wins <- floor(chrom.lens / size)  # number of windows
	warning(paste("Trimming from chromosome ends:",paste(paste(chroms,chrom.lens-size*chrom.wins,sep=": "),collapse=", "),"bp."))
	chrom.breaks <- c(0,cumsum(chrom.wins)) # 0, and then indices of *last* windows in each chromosome
    pos.fn <- function (n) {
        # return chromsome, start, and end of these windows
		this.chrom <- findInterval(n-1,chrom.breaks)
        this.chrom[ n<1 || n>max(chrom.breaks) ] <- NA
		cn <- n - chrom.breaks[this.chrom]
		win.start <- chrom.starts[this.chrom] + (cn-1)*size
		win.end <- win.start + size-1
        return( data.frame( chrom=chroms[this.chrom], start=win.start, end=win.end ) )
    }
	win.fn <- function (n,...) {
		if (n<1 || n>max(chrom.breaks)) { stop("No such window.") }
        regions <- pos.fn(n)
		vcf_query( file=file, regions=regions, samples=samples, ... )
	}
	attr(win.fn,"max.n") <- max(chrom.breaks)
	attr(win.fn,"region") <- pos.fn
    attr(win.fn,"samples") <- samples
    class(win.fn) <- c("winfun", "function")
	return(win.fn)
}

vcf_windower_snp <- function (file, sites, size, samples=vcf_samples(file)) {
	chroms <- names(sites)
	chrom.lens <- sapply( sites, length )
	chrom.wins <- floor(chrom.lens / size)
	warning(paste("Trimming from chromosome ends:",paste(paste(chroms,chrom.lens-size*chrom.wins,sep=": "),collapse=", "),"SNPs."))
	chrom.breaks <- c(0,cumsum(chrom.wins))  # 0, and then indices of *last* windows in each chromosome
    pos.fn <- function (n) {
        # return chromsome, start, and end of these windows
		this.chrom <- findInterval(n-1,chrom.breaks)
        this.chrom[ n<1 || n>max(chrom.breaks) ] <- NA
		cn <- n - chrom.breaks[this.chrom]
        win.start <- sapply( seq_along(n), function (k) { sites[[this.chrom[k]]][1+(cn[k]-1)*size] } )
        win.end <- sapply( seq_along(n), function (k) { sites[[this.chrom[k]]][cn[k]*size] } )
        return( data.frame(
                    chrom=chroms[this.chrom],
                    start=win.start,
                    end=win.end
                ) )
    }
	win.fn <- function (n,...) {
		if (n<1 || n>max(chrom.breaks)) { stop("No such window.") }
		regions <- pos.fn(n)
		vcf_query( file=file, regions=regions, samples=samples, ... )
	}
	attr(win.fn,"max.n") <- max(chrom.breaks)
	attr(win.fn,"region") <- pos.fn
    attr(win.fn,"samples") <- samples
    class(win.fn) <- c("winfun", "function")
	return(win.fn)
}

#' Window Extractor Functions
#'
#' A window extractor function ("winfun") takes an integer vector and returns an integer matrix, 
#' giving the number of alternate alleles for each sample for the corresponding windows.
#' Such functions also have three attributes: \code{max.n}, giving the index of the largest window,
#' \code{samples}, the sample IDs that are extracted (corresponding to the columns of the matrix that is returned),
#' and \code{region}, which is a function that takes an integer vector and returns a data frame giving chromosome, start, and end of the corresponding windows.
#'
#' @param f A function (for as.winfun), or a window extractor function (class \code{winfun}).
#' @param max.n Index of the largest window.
#' @param samples Character vector of sample IDs corresponding to columns of extracted data.
#' @param region A function taking an integer vector returning the chromosome, start, and end of the corresponding windows.
as.winfun <- function (f,max.n,samples,region) {
    attr(f,"max.n") <- max.n
    attr(f,"samples") <- samples
    attr(f,"region") <- region
    return(f)
}

#' @describeIn as.winfun Returns the \code{region} function of a window extractor function.
#' @export
region <- function (f) { function (n=seq_len(attr(f,"max.n"))) { attr(f,"region")(n) } }

#' @describeIn as.winfun Gives the sample IDs returned by a window extractor function.
#' @export
samples <- function (f) { attr(f,"samples") }


#' Number of Samples Returned by a Window Extractor Function
#' 
#' Gives the number of samples of matrices returned by a window extractor function (the number of rows is NA).
#'
#' @export 
#' @method dim winfun
dim.winfun <- function (f) { c(NA,length(attr(f,"samples"))) }
