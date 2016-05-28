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

#' Construct the Covariance Matrix from the Data in Pieces
#'
#' This will compute the covariance matrix for a dataset
#' that is only available in disjoint pieces (because it is too big to keep in memory at once, for instance).
#'
#' Suppose the entire data matrix is D, and n is the logical vector indicating when both D[,i] and D[,j] are nonmissing:
#'    cov(D)[i,j] = sum( D[nij,i]*D[nij,j] )/(sum(nij)-1) - sum(D[nij,i])*sum(D[nij,i])/( sum(nij)*(sum(nij)-1) )
#' @param f A function such that f(k) returns the k-th chunk of data.
#' @param n The number of chunks of data, or a vector of indices.
#' @export
running_cov <- function (f,n) {
    if (is.numeric(n) && length(n)==1) { n <- seq_len(n) }
    x <- f(1)
    colm <- colMeans(x,na.rm=TRUE)  # estimate of the column means
    x <- sweep(x,2,colm,"-")
    nn <- crossprod(!is.na(x))      # matrix of number of shared nonmissings
    mu <- colMeans(x,na.rm=TRUE)
    x[is.na(x)] <- 0
    xp <- crossprod(x)/(nn-1)       # matrix of mean product of shared nonmissings
    for (k in n) {
        x <- sweep( f(k), 2, colm, "-" )
        new.nn <- crossprod(!is.na(x))
        nn <- nn + new.nn
        mu <- ( diag(new.nn)/diag(nn) ) * mu + (1/diag(nn)) * colSums(x,na.rm=TRUE)
        x[is.na(x)] <- 0
        xp <- ( new.nn/(nn-1) ) * xp + (1/(nn-1)) * crossprod(x)/(nn-1)
    }
    return( xp - (nn/(nn-1))*outer(mu,mu,"*") )
}

