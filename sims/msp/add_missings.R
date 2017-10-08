library(lostruct)
library(parallel)

nindivs <- 1000

bcfs <- list.files(".", "*bcf$")
indiv_mut_rates <- rep(seq_len(nindivs) * 1e-6, each=2)
rel_chrom_mut_rates <- seq_along(bcfs)

mclapply(bcfs, function (bcf){
    tempfile <- paste0(bcf,"_temp_geno.txt")
    outfile <- gsub(".bcf", "_mut.bcf", bcf)
    mut_rates <- indiv_mut_rates * rel_chrom_mut_rates[match(bcf,bcfs)]
    f <- vcf_windower(file=bcf, size=1e3, type='snp')
    # write new genotypes
    for (k in seq_len(attr(f,"max.n"))) {
        x <- f(k)
        muts <- rbinom(length(x), size=1, prob=mut_rates[col(x)])
        x <- (x + muts) %% 2
        write.table(x, file=tempfile, append=TRUE, col.names=FALSE, row.names=FALSE, sep='\t')
    }
    system(paste("(bcftools view -h", bcf, 
                 "; bcftools view -H -G", bcf, "| sed -e 's/$/\tGT/' |",
                 "paste -", tempfile, ") | bcftools convert -O b - >", outfile))
    file.remove(tempfile)
  }, mc.cores=12)
