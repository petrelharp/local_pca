#' Find Eigenvectors and Eigenvalues for Covariance Matrices in Windows
#'
#' Subtracts row means from the input matrix, then
#' divides the rows into windows of a given size,
#' computes the covariance matrix for each,
#' and returns the given number of eigenvector/eigenvalue pairs.
#' Omits the last (short) window.
#'
#' @param data Either a numeric matrix or a window extractor function, that takes a positive integer \code{n} and returns the numeric matrix that correspond to the \code{n}th window.
#' @param k Number of eigenvalue/eigenvector pairs to return.
#' @param do.windows Vector of integers that correspond to the windows used.
#' @param win If \code{data} is a matrix, each contiguous block of \code{win} rows of the matrix \code{data} will be one window.
#' @param mc.cores If this is greater than 1, parallel::mclapply will be used.
#' @param ... Other parameters to be passed to \code{win.fn}.
#' @return A named list, with components
#'   $values : A numeric matrix of eigenvalues, one column for each window and k rows.
#'   $vectors : A numeric matrix, one column for each windows, having the eigenvectors in order (so it has k * ncol(matrix) rows).
#' @export
eigen_windows <- function ( 
            data, k,
            do.windows=NULL,
            win, mc.cores=1, ... ) {
    this.lapply <- if (mc.cores>1) { function (...) parallel::mclapply(...,mc.cores=mc.cores) } else { lapply }
    eigen.mat <- if (inherits(data,"function")) {
                eigen_windows_winfn( win.fn=data, do.windows=do.windows, k=k, mc.cores=mc.cores, ... )
            } else {
                if (missing(win)) { stop("Must supply either a function or a matrix and a value for win.") }
                eigen_windows_matrix( data=data, win=win, do.windows=do.windows, k=k, mc.cores=mc.cores )
            }
    return( list(
            values = eigen.mat[1:k,],
            vectors = eigen.mat[-(1:k),]
        ) )
}

#' Sets up eigen_windows_winfn to work on a matrix.
eigen_windows_matrix <- function (
            data, k, win
            do.windows=NULL,
            mc.cores ) {
    win.fn <- function(n,...){ data[seq( ((n-1)*win+1), (n*win) ),] }
    attr(win.fn,"max.n") <- floor( nrow(data)/win )
    return( eigen_windows_winfn( win.fn, do.windows=do.windows, k=k, mc.cores=mc.cores ) )
}

#' Does the work of eigen_windows.
eigen_windows_winfn <- function ( 
            win.fn=function(data,n,...){ data[seq( ((n-1)*win+1), (n*win) ),] }, 
            do.windows=NULL,
            k, mc.cores=1, ... ) {
    this.lapply <- if (mc.cores>1) { function (...) parallel::mclapply(...,mc.cores=mc.cores) } else { lapply }
    if (is.null(do.windows)) { do.windows <- seq_len(attr(win.fn,"max.n")) }
    .local <- function(x) {
        chunk <- win.fn(n=x,...)
        chunk <- sweep( chunk, 1, rowMeans(chunk,na.rm=TRUE) )
        cov <- cov(chunk,use="pairwise")
        if(any(is.na(cov))) {return(rep(NA,2*(nrow(cov)+1)))}
        PCA <- eigen(cov)
        # returns in order (values, vectors)
        return( c( PCA$values[1:k], PCA$vectors[,1:k] ) )
    }
    eigen.mat <- do.call( cbind, this.lapply( do.windows, .local ) )
    return( eigen.mat )
}
