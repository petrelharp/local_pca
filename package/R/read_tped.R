#' Read Chromosome(s) From a .tped File
#'
#' Reads specified chromsomes from a .tped file and recodes the data into numeric format:
#'    0 = major homozygote
#'    1 = heterozygote
#'    2 = minor homozygote
#' where "major" is the most common allele at that site, and all others are lumped as "minor".
#' Assumes that anything not in "ACGT" is NA (case-sensitive).
#'
#' A .tped has one row per SNP, *two* columns per sample.  The first column is chromosome, the fourth is position.
#' Does no error checking: let the user beware.
#'
#' @param file Input .tped file.
#' @param chrom Vector of chromosome name(s) to include.
#' @param triallelic Include triallelic sites? Defaults to TRUE.
#' @return An integer matrix with one row per site and one column per individual.
#' @export
read_tped <- function (file, chrom, triallelic=TRUE) {
    chr <- read.table(pipe(paste0("zcat ", file, " | grep '^\\(", paste(chrom,collapse="\\|"),"\\)\\>'")), stringsAsFactors=FALSE)
    bases <- c("A","C","G","T")
    geno <- match(as.matrix(chr[,-(1:4)]),bases)
    dim(geno) <- c( nrow(chr), ncol(chr)-4L )
    # totals[j,2] gives the total number of C's at the j-th site
    totals <- sapply( seq_along(bases), function (k) { rowSums( geno==k, na.rm=TRUE ) } )
    if (!triallelic) {
        # remove triallelic sites
        which.triallelic <- ( rowSums(totals>0) > 2 )
        totals[which.triallelic,] <- NA
    }
    maxcounts <- do.call( pmax, lapply(1:ncol(totals),function(k)totals[,k]) )
    major <- rep(NA,nrow(geno))
    for (k in seq_along(bases)) { major[ totals[,k] == maxcounts ] <- k }

    # code the genotype matrix
    coded <- matrix( NA, nrow=nrow(geno), ncol(geno) )
    coded[ geno==major ] <- 0L
    coded[ geno!=major ] <- 1L

    # convert to diploids
    n <- ncol(coded)/2
    coded <- coded[,2*(1:n)]+coded[,2*(1:n)-1]
    return(coded)
}
