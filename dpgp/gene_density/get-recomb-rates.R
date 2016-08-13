# paste into http://petrov.stanford.edu/cgi-bin/recombination-rates_updateR5.pl?results=0&release=5.36&query_type=batch&Submit=Submit
# The resulting recomb rates are in cM/Mbp.
chroms <- c("2L","2R","3L","3R","X")
lens <- structure(c(22963456L, 21142841L, 24536634L, 27894163L, 22418422L), .Dim = 5L, .Dimnames = list(c("2L", "2R", "3L", "3R", "X")))

winlen <- 1e5

for (this.chrom in chroms) {
    breaks <- c(0, winlen * (1:floor(lens[this.chrom]/winlen)), lens[this.chrom])
    writeLines(
            sprintf("%s:%d..%d",this.chrom,breaks[-length(breaks)],breaks[-1]-1),
            con=paste0("chr",this.chrom,".query")
        )
}
