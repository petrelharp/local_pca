# verify what cov(use="pairwise") does

x <- structure(
               c(2L, NA, 8L, 8L, 
                 NA, NA, 8L, 6L, 
                 NA, 1L, 9L, 2L), 
           .Dim = c(4L, 3L))

cov.x <- structure(
                   c(12,  0,  0, 
                      0,  2,  7, 
                      0,  7,  19), 
           .Dim = c(3L, 3L))

expect_equal( cov(x,use='pairwise'), cov.x )

cov_fn <- function (x,y=x) { 
    ut <- (!is.na(x)) & (!is.na(y))
    u <- x[ut]; 
    v <- y[ut]; 
    sum( (u-mean(u))*(v-mean(v)) )/(sum(ut)-1)
}

cov_fn_2 <- function (x,y=x) { 
    ut <- (!is.na(x)) & (!is.na(y))
    u <- x[ut]; 
    v <- y[ut]; 
    sum( u*v )/(sum(ut)-1) - sum(u)*sum(v)/(sum(ut)*(sum(ut)-1))
}

my.cov.2 <- my.cov <- matrix(0,nrow=nrow(cov.x),ncol=ncol(cov.x))
for (i in 1:ncol(x)) {
    for (j in 1:ncol(x)) {
        my.cov[i,j] <- cov_fn( x[,i], x[,j] )
        my.cov.2[i,j] <- cov_fn_2( x[,i], x[,j] )
    }
}

expect_equal( my.cov, cov.x )
expect_equal( my.cov.2, cov.x )
