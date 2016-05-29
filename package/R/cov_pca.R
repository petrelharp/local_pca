#' Find the Top k Eigenvalues and Eigenvectors of the Covariance Matrix
#'
#' This takes a numeric matrix, subtracts row means, finds the covariance matrix of the columns (using \code{use="pairwise"}),
#' and returns a vector containing: the sum of squared entries in the covariance matrix; the top k eigenvalues, and the top k eigenvectors.
#'
#' @param x A numeric matrix, possibly with missing values.
#' @param k Number of eigenvalue/eigenvector pairs to return.
#' @param w A vector of weights for the columns of \code{x}.
#' @param covmat The covariance matrix of \code{x}, if pre-computed (say, with \code{running_cov}).
#' @return A numeric vector: 
#' the first entry gives the total sum of squared values of the covariance matrix 
#' (i.e., the sum of the eigenvectors, for computing the proportion of variance explained);
#' the next k entires give eigenvalues and each subsequent group of ncol(data) columns give the corresponding eigenvector.
#'
#' If C is the covariance matrix, the columns of U are the eigenvectors, and W and Lambda are diagonal matrices with w and the eigenvalues on the diagonals respectively,
#' then the output is the minimizer of
#'       | W^(1/2) ( C - U Lambda U' ) W^(1/2) |^2 = sum_{ij} w[i] w[j] ( C[i,j] - (U Lambda U')[i,j] )^2 ,
#' and the sum-of-squares (equal to the sum of the eigenvalues squared) is
#'       sum_{ij}  w[i] * w[j] * C[i,j]^2 .
#' @export
cov_pca <- function (x,k,w=1,
                     covmat=cov(sweep(x,1,rowMeans(x,na.rm=TRUE),"-"),use='pairwise') ) {
    sqrt.w <- rep_len(sqrt(w),ncol(covmat))
    covmat <- sweep( sweep( covmat, 1, sqrt.w, "*" ), 2, sqrt.w, "*" )
    if(any(is.na(covmat))) {return(rep(NA,2*(nrow(covmat)+1)))}
    # PCA <- eigen(covmat)
    # return( c( sum(covmat^2), PCA$values[1:k], PCA$vectors[,1:k] ) )
    PCA <- if (k==ncol(covmat)) { eigen(covmat) } else { RSpectra::eigs_sym(covmat,k=k) }
    PCA$vectors <- sweep( PCA$vectors, 1, sqrt.w, "/" )
    # returns in order (total sumsq, values, vectors)
    return( c( sum(covmat^2), PCA$values, PCA$vectors ) )
}

#' Construct the Covariance Matrix from the Data in Pieces
#'
#' This will compute the covariance matrix for a dataset
#' that is only available in disjoint pieces (because it is too big to keep in memory at once, for instance).
#'
#' @param f A function such that f(k) returns the k-th chunk of data.
#' @param n The number of chunks of data, or a vector of indices.
#' @return A numeric matrix:
#' If a and b are vectors, then
#'    mean(a[1:(n+m)]) == (n/(n+m)) * mean(a[1:n]) + (1/(n+m)) * sum(a[(n+1):(n+m)])
#' and
#'    cov(a,b) = sum(a*b)/(n-1) - sum(a)*sum(b)/(n*(n-1))
#' where the sums are over the n shared nonmissings.
#' @export
running_cov <- function (f,n) {
    if (is.numeric(n) && length(n)==1) { n <- seq_len(n) }
    x <- f(n[1])
    z <- !is.na(x)
    colm <- colMeans(x,na.rm=TRUE)  # estimate of the column means
    x <- sweep(x,2,colm,"-")
    x[!z] <- 0
    nn <- crossprod(z)       # matrix of number of shared nonmissings
    sums <- crossprod(x,z)       # matrix of conditional sums: sums[i,j] = sum( x[,i] * !is.na(x[,j]) )
    sumsq <- crossprod(x)        # matrix of sum of products of shared nonmissings
    for (k in n[-1]) {
        x <- sweep( f(k), 2, colm, "-" )
        z <- !is.na(x)
        x[!z] <- 0
        nn <- nn + crossprod(z)
        sums <- sums + crossprod(x,z)
        sumsq <- sumsq + crossprod(x)
    }
    return( (1/(nn-1))*sumsq - sums*t(sums)/(nn*(nn-1)) )
}

