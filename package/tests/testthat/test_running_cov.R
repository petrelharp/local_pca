context("running covariance")

xx <- matrix( 1:27, ncol=3 )

f <- function (k) {
    xx[ (k-1)*3+(1:3), ]
}

n <- 1:3

expect_equal( xx, do.call(rbind, lapply(n,f) ) )


