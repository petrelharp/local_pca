context("reading windows from VCF files")
## FIX: Test.vcf has '0|0' genotypes; but query_vcf assumes '0/0'.

# has 7 individuals at 71 loci across 4334 bases
vcf.file <- "test.bcf"
bcf.file <- "test.bcf"

# read in data
#  (will warn about truncated SNP)
suppressWarnings( 
                 win.fn <- vcf_windower(bcf.file, size=7, type='snp' ) )

# get regions
these.regions <- region(win.fn)()

# do local PCA
pca.stuff <- eigen_windows( win.fn, k=4 )

