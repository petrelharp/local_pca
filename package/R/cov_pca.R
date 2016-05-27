#' Find the Top k Eigenvalues and Eigenvectors of the Covariance Matrix
#'
#' This takes a numeric matrix, subtracts row means, finds the covariance matrix of the columns (using \code{use="pairwise"}),
#' and returns a vector containing: the sum of squared entries in the covariance matrix; the top k eigenvalues, and the top k eigenvectors.
#'
#' @param x A numeric matrix, possibly with missing values.
#' @param k Number of eigenvalue/eigenvector pairs to return.
#' @param w Work in l2(sqrt(w)); see below.
#' @return A numeric vector: 
#' the first entry gives the total sum of squared values of the covariance matrix 
#' (i.e., the sum of the eigenvectors, for computing the proportion of variance explained);
#' the next k entires give eigenvalues and each subsequent group of ncol(data) columns give the corresponding eigenvector.
#'
#' If C is the covariance matrix, the columns of U are the eigenvectors, and W and Lambda are diagonal matrices with w and the eigenvalues on the diagonals respectively,
#' then the output is the least-squares solution to
#'       | W^(1/2) ( C - U Lambda U' ) W^(1/2) |,
#' and the sum-of-squares (equal to the sum of the eigenvalues squared) is
#'       sum_{ij}  ( w[i] w[j] )^(1/2) C[i,j]^2 .
#' @export
cov_pca <- function (x,k,w=1) {
    x <- sweep( x, 1, rowMeans(x,na.rm=TRUE) )
    sqrt.w <- rep_len(sqrt(w),ncol(x))
    covmat <- cov(x,use="pairwise")
    covmat <- sweep( sweep( covmat, 1, sqrt.w, "*" ), 2, sqrt.w, "*" )
    if(any(is.na(covmat))) {return(rep(NA,2*(nrow(covmat)+1)))}
    # PCA <- eigen(covmat)
    # return( c( sum(covmat^2), PCA$values[1:k], PCA$vectors[,1:k] ) )
    PCA <- if (k==ncol(x)) { eigen(covmat) } else { RSpectra::eigs_sym(covmat,k=k) }
    PCA$vectors <- sweep( PCA$vectors, 1, sqrt.w, "/" )
    # returns in order (total sumsq, values, vectors)
    return( c( sum(covmat^2), PCA$values, PCA$vectors ) )
}

