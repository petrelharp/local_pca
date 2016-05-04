#' Find Eigenvectors and Eigenvalues for Covariance Matrices in Windows
#'
#' Subtracts row means from the input matrix, then
#' divides the rows into windows of a given size,
#' computes the covariance matrix for each,
#' and returns the given number of eigenvector/eigenvalue pairs.
#' Omits the last (short) window.
#'
#' @param data The object containing the genotype data, e.g., a numeric matrix.
#' @param win.fn Function that takes \code{data} and a positive integer \code{n} and returns the numeric matrix that correspond to the \code{n}th window.
#' @param do.windows Vector of integers (or, arguments to \code{win.fn}) that correspond to the windows used.
#' @param win If \code{win.fn} is not provided, each contiguous block of \code{win} rows of the matrix \code{data} will be one window.
#' @param k Number of eigenvalue/eigenvector pairs to return.
#' @param mc.cores If this is greater than 1, parallel::mclapply will be used.
#' @param ... Other parameters to be passed to \code{win.fn}.
#' @return A named list, with components
#'   $values : A numeric matrix of eigenvalues, one column for each window and k rows.
#'   $vectors : A numeric matrix, one column for each windows, having the eigenvectors in order (so it has k * ncol(matrix) rows).
#' @export
eigen_windows <- function ( 
            data, 
            win.fn=function(data,n,...){ data[seq( ((n-1)*win+1), (n*win) ),] }, 
            do.windows=1:floor(nrow(data)/win),
            win, k, mc.cores=1, ... ) {
    if (missing(win) && missing(win.fn)) { stop("Must supply either win or win.fn.") }
    this.lapply <- if (mc.cores>1) { function (...) parallel::mclapply(...,mc.cores=mc.cores) } else { lapply }
    .local <- function(x) {
        chunk <- win.fn(data,n=x,...)
        chunk <- sweep( chunk, 1, rowMeans(chunk,na.rm=TRUE) )
        cov <- cov(chunk,use="pairwise")
        if(any(is.na(cov))) {return(rep(NA,2*(nrow(cov)+1)))}
        PCA <- eigen(cov)
        # returns in order (values, vectors)
        return( c( PCA$values[1:k], PCA$vectors[,1:k] ) )
    }
    eigen.mat <- do.call( cbind, this.lapply( do.windows, .local ) )
    return( list(
            values = eigen.mat[1:k,],
            vectors = eigen.mat[-(1:k),]
        ) )
}
