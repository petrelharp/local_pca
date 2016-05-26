context("reading windows from VCF files")


# has 7 individuals at 71 loci across 4334 bases
bcf.file <- "test.vcf"

# read in data
win.fn <- vcf_windower(bcf.file, size=7, type='snp' )

# get regions
these.regions <- region(win.fn)()

# do local PCA
pca.stuff <- eigen_windows( win.fn, k=4 )

