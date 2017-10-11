context("running covariance")

set.seed(123)
nchunks <- 100
chunklen <- 5
xx <- matrix( rpois(nchunks*chunklen*3, 5), ncol=3 )

f <- function (k) {
    xx[ (k-1)*chunklen+(1:chunklen), ]
}

n <- 1:nchunks

expect_equal( xx, do.call(rbind, lapply(n,f) ) )

expect_equal( cov(xx,use="pairwise"), running_cov(f,n, normalize.rows=FALSE) )

## with normalize.rows

expect_equal( cov(sweep(xx,1,rowMeans(xx)),use="pairwise"), running_cov(f,n, normalize.rows=TRUE) )

Imat <- diag(3) - 1/3
expect_equal( Imat%*%cov(xx,use="pairwise")%*%Imat, running_cov(f,n, normalize.rows=TRUE) )

## with missings

xx[1,2] <- xx[2,2] <- xx[1,3] <- xx[6,1] <- xx[7,2] <- NA

expect_equal( xx, do.call(rbind, lapply(n,f) ) )

expect_equal( cov(xx,use="pairwise"), running_cov(f,n, normalize.rows=FALSE) )
