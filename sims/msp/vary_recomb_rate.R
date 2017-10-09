library(lostruct)
library(parallel)

bcfs <- list.files(".", "*bcf$")
rel_genetic_dist <- seq(1, 10, length.out=length(bcfs))

mclapply(bcfs, function (bcf){
    tempfile <- paste0(bcf,"_temp_map.txt")
    outfile <- gsub(".bcf", "_recomb.bcf", bcf)
    f <- vcf_windower(file=bcf, size=1e3, type='snp')
    system(paste("(bcftools view -h", bcf, 
                 "; bcftools view -H", bcf, 
                 "| awk 'BEGIN {OFS=\"\\t\"} {$2 = sprintf(\"%.0f\", $2*2.5)}1') | bcftools convert -O b - >", outfile))
    file.remove(tempfile)
  }, mc.cores=12)

