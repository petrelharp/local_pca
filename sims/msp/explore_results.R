library(lostruct)
library(parallel)

basedir <- "nada_015374/flat_mut_024499"

bcf_files <- list.files(basedir, ".*bcf$", full.names=TRUE)

winfns <- lapply(bcf_files, vcf_windower, size=1e6, type='bp')
