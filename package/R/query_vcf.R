#' Find Minimal Enclosing Regions
#'
#' Given a named list of positions, whose names correspond to chromosomes,
#' return the data frame of regions with names chrom,start,end that encloses them.
#'
#' @param sites Named list of positions, as output by \code{vcf_positions()}.
#' @return A data frame with variables chrom, start, and end.
#' @export
enclosing_region <- function (sites) {
    if (length(sites)==0 || max(sapply(sites,length))==0) { return( NULL ) }
    return( data.frame(
            chrom=names(sites),
            start=sapply(sites,min,na.rm=TRUE),
            end=sapply(sites,max,na.rm=TRUE)
        ) )
}

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
#' returning the genotypes in numeric format.
#' Assumes the data are diploid.
#'
#' @param file The name of the indexed vcf file.
#' @param regions A data frame containing chrom, start, and end of the regions to be extracted.
#' @param samples A character vector of sample names to be extracted.
#' @param verbose Whether to output with \code{cat()} the bcftools calls.
#' @return An integer matrix, giving the number of alternate alleles for each sample.
#' @export
query_genotypes <- function (file,regions,samples, verbose=FALSE) {
	bcf.args <- c("bcftools", "query", "-f", "'[ %GT]\\n'", "-r", region_string(regions))
	if (!missing(samples) && length(samples)>0) { bcf.args <- c( bcf.args, "-s", paste(samples,collapse=',') ) }
	bcf.args <- c( bcf.args, file )
    bcf.call <- paste(bcf.args, collapse=" ")
    if (verbose) cat(bcf.call, "\n")
    gt.text <- data.table::fread( bcf.call, header=FALSE, sep=' ', data.table=FALSE )
    gt <- c(0L,1L,1L,2L)[match( unlist(gt.text), c("0/0","0/1","1/0","1/1") )]
    dim(gt) <- dim(gt.text)
    return(gt)
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

#' Construct a Window Extractor for a VCF File
#'
#' Returns a function that, given an integer \code{n}, returns a numeric matrix of genotypes from the the \code{n}th window of the VCF file.
#'
#' @param file The name of the indexed vcf file.
#' @param size The size of the window.
#' @param type The units of the window: 'bp' or 'snp'?
#' @param sites The positions in the VCF file, as returned by \code{vcf_positions}.
#' @param samples A character vector of sample names to be extracted.
#' @param f A window extractor function.
#' @return A class "winfun" function that returns an integer matrix, giving the number of alternate alleles for each sample.
#' Such functions also have two attributes: \code{max.n}, giving the index of the largest window,
#' and \code{region}, which is a function that takes an integer vector and returns a data frame giving chromosome, start, and end of the corresponding windows.
#' @export
vcf_windower <- function (
			file, 
			size, 
			type, 
			sites=vcf_positions(file),
			samples=NULL) {
	if (type=="bp") {
		vcf_windower_bp( file=file, sites=sites, size=size, samples=samples )
	} else if (type=="snp") {
		vcf_windower_snp( file=file, sites=sites, size=size, samples=samples )
	}
}

vcf_windower_bp <- function (file, sites, size, samples=NULL) {
	chroms <- names(sites)
	chrom.lens <- sapply( sites, max )
	chrom.wins <- floor(chrom.lens / size)
	warning(paste("Trimming from chromosome ends:",paste(paste(chroms,chrom.lens-size*chrom.wins,sep=": "),collapse=", "),"bp."))
	chrom.breaks <- c(0,cumsum(chrom.wins))
    pos.fn <- function (n) {
        # return chromsome, start, and end of these windows
		this.chrom <- findInterval(n-1,chrom.breaks)
        this.chrom[ n<1 || n>max(chrom.breaks) ] <- NA
		cn <- n - chrom.breaks[this.chrom]
		win.start <- 1+(cn-1)*size
		win.end <- win.start + size-1
        return( data.frame( chrom=chroms[this.chrom], start=win.start, end=win.end ) )
    }
	win.fn <- function (n) {
		if (n<1 || n>max(chrom.breaks)) { stop("No such window.") }
        regions <- pos.fn(n)
		query_genotypes( file=file, regions=regions, samples=samples )
	}
	attr(win.fn,"max.n") <- max(chrom.breaks)
	attr(win.fn,"region") <- pos.fn
    class(win.fn) <- c("winfun", "function")
	return(win.fn)
}

vcf_windower_snp <- function (file, sites, size, samples=NULL) {
	chroms <- names(sites)
	chrom.lens <- sapply( sites, length )
	chrom.wins <- floor(chrom.lens / size)
	warning(paste("Trimming from chromosome ends:",paste(paste(chroms,chrom.lens-size*chrom.wins,sep=": "),collapse=", "),"SNPs."))
	chrom.breaks <- c(0,cumsum(chrom.wins))  # these are indices of *last* windows in each chromosome
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
	win.fn <- function (n) {
		if (n<1 || n>max(chrom.breaks)) { stop("No such window.") }
		regions <- pos.fn(n)
		query_genotypes( file=file, regions=regions, samples=samples )
	}
	attr(win.fn,"max.n") <- max(chrom.breaks)
	attr(win.fn,"region") <- pos.fn
    class(win.fn) <- c("winfun", "function")
	return(win.fn)
}

#' @describeIn vcf_windower Returns the \code{region} function of a window extractor function.
region <- function (f) {
    function (n=seq_len(attr(f,"max.n"))) { attr(f,"region")(n) }
}
