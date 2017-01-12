context("reading windows from VCF files")
## FIX: Test.vcf has '0|0' genotypes; but query_vcf assumes '0/0'.

# has 7 individuals at 71 loci across 4334 bases
# need to run
#   bcftools view -O b -o test.bcf test.vcf
#   bcftools index test.bcf
#  
bcf.file <- "test.bcf"

vcf.text <- read.table("test.vcf", sep="\t", skip=30, header=TRUE, comment.char="" )
vcf.mat <- vcf.text[,10:16]

# read in data
#  (will warn about truncated SNP)
expect_warning( 
                 win.fn <- vcf_windower(bcf.file, size=7, type='snp' ) 
         )
expect_equal( attr(win.fn,"max.n"), 10 )

# get regions
expect_equal( region(win.fn)(),
             structure(list(chrom = structure(c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                                1L, 1L, 1L), .Label = "2", class = "factor"), start = c(10179L, 
                                              10199L, 10574L, 13630L, 13811L, 13970L, 14058L, 14151L, 14238L, 
                                              14411L), end = c(10190L, 10572L, 10753L, 13750L, 13949L, 14027L, 
                                              14146L, 14234L, 14393L, 14494L)), .Names = c("chrom", "start", 
                            "end"), row.names = c(NA, -10L), class = "data.frame")
             )

expect_equal( colnames(vcf.mat), samples(win.fn) )

# get the first window
expect_equal( win.fn(1,recode=FALSE),
            structure(c("0|0", "0|0", "0|0", "0|0", "0|0", "0|0", "0|0", 
            "0|0", "0|0", "0|0", "0|0", "0|1", "0|0", "0|0", "0|0", "0|0", 
            "0|0", "0|0", "2|0", "0|0", ".", "0|0", "0|0", "0|0", "0|0", 
            "0|0", "0|0", "0|0", "0|0", "0|0", "0|0", "0|0", "1|2", "0|0", 
            "1|1", "0|0", "0|0", "0|1", "0|0", "0|0", "0|0", "1|0", "0|0", 
            "0|0", "0|0", "0|0", "0|0", "0|0", "0|0"), .Dim = c(7L, 7L), .Dimnames = list(
                NULL, c("V1", "V2", "V3", "V4", "V5", "V6", "V7")))
        )

expect_equivalent( win.fn(1,recode=FALSE), as.matrix(vcf.mat[1:7,]) )


expect_equal( win.fn(1), 
            structure(c(0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 
            0L, 0L, 0L, 0L, 0L, 1L, 0L, NA, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
            0L, 0L, 0L, 2L, 0L, 2L, 0L, 0L, 1L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 
            0L, 0L, 0L, 0L), .Dim = c(7L, 7L))
    )

expect_equal( win.fn(2),
            structure(c(0L, 0L, NA, 0L, 0L, 0L, 0L, 0L, 0L, NA, 0L, 0L, 0L, 
            0L, NA, 0L, NA, 0L, 1L, 0L, 0L, NA, 0L, NA, 1L, 1L, 0L, 0L, 0L, 
            0L, NA, 0L, 0L, 0L, 1L, 0L, 0L, NA, 0L, 0L, 0L, 2L, 0L, 0L, NA, 
            0L, 0L, 0L, 0L), .Dim = c(7L, 7L))
    )


# do local PCA
pca.stuff <- eigen_windows( win.fn, k=2 )

expect_equivalent( pca.stuff[1,2:3],  eigen( cov(sweep( win.fn(1), 1, rowMeans(win.fn(1),na.rm=TRUE), "-" ), use='pairwise') )$values[1:2] )
expect_equivalent( pca.stuff[2,2:3],  eigen( cov(sweep( win.fn(2), 1, rowMeans(win.fn(2),na.rm=TRUE), "-" ), use='pairwise') )$values[1:2] )

# now with bp 

expect_warning( 
                 win.fn <- vcf_windower(bcf.file, size=400, type='bp' ) 
         )

expect_equal( attr(win.fn,"max.n"), 10 )

# get regions
expect_equal( region(win.fn)(),
            structure(list(chrom = structure(c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 
            1L, 1L, 1L), .Label = "2", class = "factor"), start = c(10179, 
            10579, 10979, 11379, 11779, 12179, 12579, 12979, 13379, 13779
            ), end = c(10578, 10978, 11378, 11778, 12178, 12578, 12978, 13378, 
            13778, 14178)), .Names = c("chrom", "start", "end"), row.names = c(NA, 
            -10L), class = "data.frame")
    )

expect_true( all( region(win.fn)()[,"end"] - region(win.fn)()[,"start"] == 399 ) )

expect_equal( 
        sapply( lapply( 1:10, win.fn ), length ),
        c(105L, 42L, 0L, 0L, 0L, 0L, 0L, 0L, 49L, 154L)
    )


expect_equal( win.fn(1), 
            structure(c(0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, NA, 0L, 0L, 0L, 
            0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, NA, 0L, 0L, 0L, 0L, 
            2L, 0L, 0L, 0L, 0L, 1L, 0L, NA, NA, 0L, NA, 0L, 1L, 0L, 0L, 0L, 
            0L, 0L, 0L, 0L, 0L, 0L, 0L, NA, 0L, NA, 1L, 1L, 0L, 0L, 0L, 0L, 
            0L, 0L, 0L, 2L, 0L, 2L, 0L, 0L, NA, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 
            1L, 0L, 0L, 0L, 1L, 0L, 0L, NA, 0L, 0L, 0L, 2L, 0L, 0L, 0L, 0L, 
            0L, 0L, 0L, 0L, 0L, 0L, NA, 0L, 0L, 0L, 0L, 0L), .Dim = c(15L, 
            7L))
    )

expect_equal( win.fn(3), NULL )


# do local PCA
pca.stuff <- eigen_windows( win.fn, k=2 )
# not enough nonzero eigenvalues
expect_warning( eigen_windows( win.fn, k=4 ) )

expect_equivalent( pca.stuff[1,2:3],  eigen( cov(sweep( win.fn(1), 1, rowMeans(win.fn(1),na.rm=TRUE), "-" ), use='pairwise') )$values[1:2] )
expect_equivalent( pca.stuff[2,2],  eigen( cov(sweep( win.fn(2), 1, rowMeans(win.fn(2),na.rm=TRUE), "-" ), use='pairwise') )$values[1] )
expect_true( all( is.na( pca.stuff[3,] ) ) )

