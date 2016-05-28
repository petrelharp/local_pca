context("running covariance")

xx <- matrix( 1:27, ncol=3 )

f <- function (k) {
    xx[ (k-1)*3+(1:3), ]
}

n <- 1:3

expect_equal( xx, do.call(rbind, lapply(n,f) ) )

expect_equal( cov(xx,use="pairwise"), running_cov(f,n) )

## with missings

xx[1,2] <- xx[2,2] <- xx[1,3] <- xx[6,1] <- xx[5,2] <- NA

expect_equal( cov(xx,use="pairwise"), running_cov(f,n) )

