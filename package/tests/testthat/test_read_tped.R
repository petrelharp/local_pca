ans <- structure(c(0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 2L, 0L, 0L, 1L, 
                   1L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 
                   1L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 
                   0L, 1L, 0L, 0L, 0L, NA, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, NA), .Dim = c(20L, 3L))

tped.mat <- read_tped("test.tped",22)
tped.gz.mat <- read_tped("test.tped.gz",22)

expect_equal( tped.mat, ans )
expect_equal( tped.gz.mat, ans )

