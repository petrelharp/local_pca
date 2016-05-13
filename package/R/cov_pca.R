#' Find the Top k Eigenvalues and Eigenvectors of the Covariance Matrix
#'
#' This takes a numeric matrix, subtracts row means, finds the covariance matrix of the columns (using \code{use="pairwise"}),
#' and returns a vector containing: the sum of squared entries in the covariance matrix; the top k eigenvalues, and the top k eigenvectors.
#'
#' @param x A numeric matrix, possibly with missing values.
#' @param k Number of eigenvalue/eigenvector pairs to return.
#' @return A numeric vector: 
#' the first entry gives the total sum of squared values of the covariance matrix 
#' (i.e., the sum of the eigenvectors, for computing the proportion of variance explained);
#' the next k entires give eigenvalues and each subsequent group of ncol(data) columns give the corresponding eigenvector.
#' @export
cov_pca <- function (x,k) {
    x <- sweep( x, 1, rowMeans(x,na.rm=TRUE) )
    covmat <- cov(x,use="pairwise")
    if(any(is.na(covmat))) {return(rep(NA,2*(nrow(covmat)+1)))}
    PCA <- eigen(covmat)
    # returns in order (values, vectors)
    return( c( sum(covmat^2), PCA$values[1:k], PCA$vectors[,1:k] ) )
}

