context("testing covariance/pca computations")

x <- structure(c(7L, 0L, 0L, 1L, 5L, 2L, 8L, 4L, 1L, 5L, 9L, 8L), .Dim = c(4L, 3L))

covmat <- cov(x-rowMeans(x,na.rm=TRUE))

cov.pca <- cov_pca(x,k=2)

# check unweighted PCA
expect_equal( cov.pca[1:3],
             c( sum(covmat^2),
                eigen(covmat)$values[1:2] ) )
for (k in 1:2) {
    expect_equal( rep(0,ncol(x)), 
                    pmin( abs(cov.pca[3+(k-1)*ncol(x)+(1:ncol(x))]-eigen(covmat)$vectors[,k]),
                          abs(cov.pca[3+(k-1)*ncol(x)+(1:ncol(x))]+eigen(covmat)$vectors[,k]) ) )
}

# reconstruction from eigenvalues
cov.pca <- cov_pca(x,k=ncol(x))
expect_equal( covmat,
             { U <- matrix( cov.pca[-(1:(1+ncol(x)))], ncol=ncol(x) ); 
               L <- diag( cov.pca[2:(1+ncol(x))] );
               U %*% L %*% t(U) } )

# and now weighted
w <- 1:ncol(x)

covmat <- cov(x-rowMeans(x,na.rm=TRUE))

cov.pca <- cov_pca(x,k=2,w=w)

expect_equal( cov.pca[1:3],
             c( sum( (outer(sqrt(w),sqrt(w),"*") * covmat)^2 ),
                eigen( outer(sqrt(w),sqrt(w),"*") * covmat )$values[1:2] ) )

# reconstruction, again
cov.pca <- cov_pca(x,k=ncol(x),w=w)
expect_equal( covmat,
             { U <- matrix( cov.pca[-(1:(1+ncol(x)))], ncol=ncol(x) ); 
               L <- diag( cov.pca[2:(1+ncol(x))] );
               U %*% L %*% t(U) } )
