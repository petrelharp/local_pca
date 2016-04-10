context("Testing reading TPED file")

ans <- structure(c(0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 2L, 0L, 0L, 1L, 
                   1L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 
                   1L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 
                   0L, 1L, 0L, 0L, 0L, NA, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, NA), .Dim = c(20L, 3L))

ans.phased <- structure(c(0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 1L, 0L, 0L, 1L, 
                1L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
                1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
                0L, 0L, 1L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
                0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
                0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 1L, 0L, 0L, 0L, NA, 0L, 0L, 
                0L, 0L, 0L, 0L, 0L, 0L, NA, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
                0L, NA, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, NA), .Dim = c(20L, 6L))

tped.mat <- read_tped("test.tped",22)
tped.gz.mat <- read_tped("test.tped.gz",22)
tped.mat.2 <- read_tped("test.tped")
tped.phased.mat <- read_tped("test.tped",22,phased=TRUE)

expect_equal( tped.mat, ans )
expect_equal( tped.mat.2, ans )
expect_equal( tped.gz.mat, ans )
expect_equal( tped.phased.mat[,c(1,3,5)] + tped.phased.mat[,c(2,4,6)], tped.mat )
expect_equal( tped.phased.mat, ans.phased )

context("Testing reading VCF file")

ans.vcf <- structure(c(NA, 0L, 0L, 0L, 0L, 0L, 0L, NA, 0L, 0L, 0L, 0L, 0L, 
                NA, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
                0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
                0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
                0L, 0L, 0L, NA, 0L, 0L, 0L, 0L, 0L, 0L, NA, 0L, 0L, 0L, 1L, 0L, 
                0L, NA, 0L, 0L, 0L, 0L, 0L, NA, 2L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
                1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
                0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
                0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, NA, 0L, 0L, 0L, 0L, 0L, 
                0L, NA, 0L, 0L, 0L, 1L, 0L, 0L, NA, 0L, 0L, 0L, 1L, 0L, NA, 0L, 
                0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
                0L, 0L, 0L, NA, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
                0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
                0L, NA, 0L, 0L, 0L, 0L, 0L, 0L, NA, 0L, 0L, 0L, 0L, 0L, 0L, NA, 
                0L, 0L, 1L, 1L, 0L, NA, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 
                0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, NA, 0L, 0L, 0L, 0L, 0L, 
                0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
                0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, NA, 0L, 0L, 0L, 0L, 0L, 0L, NA, 
                0L, 0L, 0L, 2L, 0L, NA, NA, 0L, 0L, 0L, 0L, 0L, NA, 0L, 0L, 0L, 
                0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
                0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
                0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, NA, 
                0L, 0L, 0L, 0L, 0L, 0L, NA, 0L, 1L, 0L, 0L, 0L, NA, NA, 0L, 0L, 
                0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 
                0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, NA, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
                0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
                0L, 0L, 0L, 0L, 0L, 0L, NA, 0L, 0L, 0L, 0L, 0L, 0L, NA, 0L, 0L, 
                0L, 0L, 0L, 0L, NA, 0L, 0L, 0L, 0L, 0L, NA, 0L, 0L, 0L, 0L, 0L, 
                0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, NA, 
                0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
                0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, NA, 0L, 0L, 
                0L, 0L, 0L, 0L), .Dim = c(71L, 7L))


true.vcf <- read.table("test.vcf.tsv", header=FALSE, skip=1, stringsAsFactors=FALSE)
colnames(true.vcf) <- c("REF","ALT",paste0(rep(scan("test.vcf.tsv",what="char",nlines=1)[-(1:2)],2),c("_1","_2")))
true.vcf.mat <- as.matrix(true.vcf[,-(1:2)])
true.geno <- t( sapply( 1:nrow(true.vcf.mat), function (k) {
            alleles <- c( true.vcf[k,"REF"], strsplit( true.vcf[k,"ALT"], "," )[[1]] )
            alleles[ 1+as.vector(true.vcf.mat[k,]) ]
        } ) )
true.mat <- recode_numeric(true.geno)
true.phased.mat <- recode_numeric(true.geno,ploidy=1)

vcf.mat <- read_vcf("test.vcf")
vcf.phased.mat <- read_vcf("test.vcf",phased=TRUE)

expect_equal( vcf.mat, true.mat )
expect_equal( vcf.mat, ans.vcf )
expect_equal( vcf.phased.mat, true.phased.mat )
