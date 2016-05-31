#' Find Eigenvectors and Eigenvalues for Covariance Matrices in Windows
#'
#' Subtracts row means from the input matrix, then
#' divides the rows into windows of a given size,
#' computes the covariance matrix for each,
#' and returns for each window the total sum of squares
#' and the given number of eigenvector/eigenvalue pairs.
#' Omits the last (short) window.
#'
#' @param data Either a numeric matrix or a window extractor function, that takes a positive integer \code{n} and returns the numeric matrix that correspond to the \code{n}th window.
#' @param k Number of eigenvalue/eigenvector pairs to return.
#' @param do.windows Vector of integers that correspond to the windows used.
#' @param win If \code{data} is a matrix, each contiguous block of \code{win} rows of the matrix \code{data} will be one window.
#' @param w A vector of weights corresponding to the columns of the data.
#' @param mc.cores If this is greater than 1, parallel::mclapply will be used.
#' @param ... Other parameters to be passed to \code{win.fn}.
#' @return A numeric matrix with one row for each window; the first column is the sum of squared values of the covariance matrix;
#' the next k columns are the eigenvalues, and each subsequent group of ncol(data) columns give the corresponding eigenvector.
#' The result has an attribute, "npc", which is equal to the number of eigenvalues computed (k).
#'
#' The work of finding eigenvalues/vectors is done by \code{cov_pca()}; see there for details.
#' @export
eigen_windows <- function (
            data, k,
            do.windows=NULL,
            win, 
            w=1,
            mc.cores=1, 
            ... ) {
    eigen.mat <- if (inherits(data,"function")) {
                eigen_windows_winfn( win.fn=data, do.windows=do.windows, k=k, w=w, mc.cores=mc.cores, ... )
            } else {
                if (missing(win)) { stop("Must supply either a function or a matrix and a value for win.") }
                eigen_windows_matrix( data=data, win=win, do.windows=do.windows, k=k, w=w, mc.cores=mc.cores )
            }
    return( eigen.mat )
}

#' Sets up a winfun-like accessor to a matrix.
winfun_matrix <- function (data,win) {
    win.fn <- function(n,...){ data[seq( ((n-1)*win+1), (n*win) ),] }
    attr(win.fn,"max.n") <- floor( nrow(data)/win )
    attr(win.fn,"samples") <- if (is.null(colnames(data))) { 1:ncol(data) } else { colnames(data) }
    attr(win.fn,"region") <- function (n) { data.frame( chrom="matrix", start=((n-1)*win+1), end=n*win ) }
    return( win.fn )
}

#' Sets up eigen_windows_winfn to work on a matrix.
eigen_windows_matrix <- function (
            data, k, win, w,
            do.windows=NULL,
            mc.cores ) {
    return( eigen_windows_winfn( winfun_matrix(data,win=win), do.windows=do.windows, k=k, w=w, mc.cores=mc.cores ) )
}

#' Does the work of eigen_windows.
eigen_windows_winfn <- function ( 
            win.fn,
            do.windows=NULL,
            k, w, mc.cores=1, ... ) {
    this.lapply <- if (mc.cores>1) { function (...) parallel::mclapply(...,mc.cores=mc.cores) } else { lapply }
    if (is.null(do.windows)) { do.windows <- seq_len(attr(win.fn,"max.n")) }
    samples <- samples(win.fn)
    out.colnames <- c( "total", paste0("lam_",1:k), as.vector( paste0("PC_", t(outer( 1:k, samples, paste, sep="_" )) ) ) )
    empty.value <- rep( NA, length(out.colnames) )
    .local <- function(n) {
        # returns in order (total, values, vectors)
        x <- win.fn(n=n,...)
        if (is.null(x)) {  # if there are no markers in this window
            empty.value
        } else {
            cov_pca( x, k=k, w=w )
        }
    }
    eigen.mat <- do.call( rbind, this.lapply( do.windows, .local ) )
    colnames(eigen.mat) <- out.colnames
    attr(eigen.mat,"npc") <- k
    attr(eigen.mat,"w") <- w
    return( eigen.mat )
}
