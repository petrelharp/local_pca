samples <- read.table("sample_info.tsv",sep="\t",header=TRUE)

samples$nsamples <- table(samples$Country.of.Origin)[samples$Country.of.Origin]
samples$weight <- 1/pmax(10,samples$nsamples)

bcf.samples <- vcf_samples("chr1-filtered-set-2014Apr15.bcf")
bcf.samples[bcf.samples=="HM020-I"] <- "HM020"

weights <- data.frame( ID=bcf.samples, weight=samples$weight[match(bcf.samples,samples$ID)] )

write.table(weights, sep="\t", file="inverse-samplesize-weights.tsv", row.names=FALSE)
