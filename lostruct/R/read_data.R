#' Convert Matrix of Alleles to Numeric
#'
#' Recodes a character matrix, where each row is a site and each \code{ploidy} columns are the alleles of an individual, into numeric format,
#' counting the number of major alleles, where "major" is the most common allele at that site, and all others are lumped as "minor".
#'
#' For speed, currently only deals with alleles that are A, C, G, or T (not indels).
#'
#' Values range from 0 (all minor alleles) to \code{ploidy} (all major alleles).
#'
#' @param x Character matrix of alleles.
#' @param ploidy Number of columns per individual.
#' @param triallelic Include triallelic sites? Defaults to TRUE; otherwise sets these to NA.
#' @param alleles A character vector of valid alleles.
#' @return An integer matrix with one row per site and one column per individual.
#' @export
recode_numeric <- function (x, ploidy=2, triallelic=TRUE, alleles=c("A","C","G","T")) {
    geno <- match(as.matrix(x),alleles)
    dim(geno) <- c( nrow(x), ncol(x) )
    # totals[j,2] gives the total number of C's at the j-th site
    totals <- sapply( seq_along(alleles), function (k) { rowSums( geno==k, na.rm=TRUE ) } )
    if (!triallelic) {
        # remove triallelic sites
        which.triallelic <- ( rowSums(totals>0) > 2 )
        totals[which.triallelic,] <- NA
    }
    maxcounts <- do.call( pmax, lapply(1:ncol(totals),function(k)totals[,k]) )
    major <- rep(NA,nrow(geno))
    for (k in seq_along(alleles)) { major[ totals[,k] == maxcounts ] <- k }

    # code the genotype matrix
    coded <- matrix( NA, nrow=nrow(geno), ncol(geno) )
    coded[ geno==major ] <- 0L
    coded[ geno!=major ] <- 1L

    # add up alleles by individual
    n <- ncol(coded)/ploidy
    out <- coded[,1+ploidy*(0:(n-1))]
    for (k in seq_len(ploidy)[-1]) {
        out <- out + coded[,k+ploidy*(0:(n-1))]
    }
    return(out)
}


#' Read Chromosome(s) From a TPED File
#'
#' Reads specified chromsomes from a .tped file and recodes the data into numeric format, e.g. for diploids:
#'    0 = major homozygote
#'    1 = heterozygote
#'    2 = minor homozygote
#' where "major" is the most common allele at that site, and all others are lumped as "minor".
#' Assumes that anything not in "ACGT" is NA (case-sensitive).
#'
#' A .tped has one row per SNP, *two* columns per sample.  The first column is chromosome, the fourth is position.
#' Does no error checking: let the user beware.
#'
#' @param file Input .tped filename.
#' @param chrom Vector of chromosome name(s) to include (defaults to all).
#' @param triallelic Include triallelic sites? Defaults to TRUE; otherwise sets these to NA.
#' @param phased Are the data phased? (i.e., should return one column per individual or two?)
#' @return An integer matrix with one row per site and one (or two, if phased) column per individual.
#' @export
read_tped <- function (file, chrom="[^ \t]*", triallelic=TRUE, phased=FALSE) {
    thecat <- if ( grepl(".gz$",file) ) { "zcat" } else { "cat" }
    # chr <- read.table(pipe(paste0(thecat, " ", file, " | grep '^\\(", paste(chrom,collapse="\\|"),"\\)\\>'")), stringsAsFactors=FALSE)
    chr <- data.table::fread((paste0(thecat, " ", file, " | grep '^\\(", paste(chrom,collapse="\\|"),"\\)\\>'")), stringsAsFactors=FALSE)
    return( recode_numeric( chr[,-(1:4)], ploidy=2-phased, triallelic=triallelic ) )
}

#' Read Chromosome(s) From a VCF File
#'
#' Reads specified chromsomes from a .vcf file and recodes the data into numeric format:
#'    0 = major homozygote
#'    1 = heterozygote
#'    2 = minor homozygote
#' where "major" is the most common allele at that site, and all others are lumped as "minor".
#' Assumes that anything not in "ACGT" is NA (case-sensitive).
#'
#' @param file Input .vcf file, possibly gzip'ped.
#' @param phased Are the data phased? (i.e., should return one column per individual or two?)
#' @param triallelic Include triallelic sites? Defaults to TRUE.
#' @return An integer matrix with one row per site and one column per individual.
#' @export
read_vcf <- function (file, phased=FALSE, triallelic=TRUE) {
    # # bare-bones style if all entries are GT format
    # vcf.header <- scan( pipe(paste0("zcat ", file, " | grep -v '^##' | head -n 1")), what='char')
    # thecat <- if (grepl(".gz$",file)) { "zcat" } else { "cat" }
    # chr <- data.table::fread(paste0(thecat, " ", file, " | grep -v '^#' | cut -f 9- | sed -e 's_0[|/]0_0_g' -e 's_0[|/][1-9]_1_g' -e 's_[1-9][|/]0_1_g' -e 's_[1-9][|/][1-9]_2_g' "), sep='\t', stringsAsFactors=FALSE, ...)
    ## etcetera
    dips <- vcf_query( file, recode=FALSE )
    haps <- unlist(strsplit(unlist(dips),"[/|]"))
    dim(haps) <- c(2,dim(dips))
    haps <- aperm( haps, c(2,1,3) )
    dim(haps) <- dim(dips)*c(1,2)
    return( recode_numeric( haps, ploidy=2-phased, triallelic=triallelic, alleles=as.character(0:3) ) )
}

