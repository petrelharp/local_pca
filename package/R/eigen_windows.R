#' Find Eigenvectors and Eigenvalues for Covariance Matrices in Windows
#'
#' Subtracts row means from the input matrix, then
#' divides the rows into windows of a given size,
#' computes the covariance matrix for each,
#' and returns the given number of eigenvector/eigenvalue pairs.
#' Omits the last (short) window.
#'
#' @param data A numeric matrix.
#' @param win Number of rows to take as one window.
#' @param k Number of eigenvalue/eigenvector pairs to return.
#' @param do.parallel Use mclapply?
#' @return A named list, with components
#'   $values : A numeric matrix of eigenvalues, one column for each window and k rows.
#'   $vectors : A list whose j-th element is a numeric matrix of j-th eigenvectors, one column for each window.
#' @export
eigen_windows <- function (data, win, k, do.parallel=TRUE) {
    this.lapply <- if (do.parallel) { function (...) parallel::mclapply(...,mc.cores=parallel::detectCores()) } else { lapply }
    data <- sweep( data, 1, rowMeans(data,na.rm=TRUE) )
    .local <- function(x) {
        chunk <- data[((x-1)*win + 1):(x*win), ]
        cov <- cov(chunk,use="pairwise")
        if(any(is.na(cov))) {return(rep(NA,2*(nrow(cov)+1)))}
        PCA <- eigen(cov)
        # returns in order (values, vectors)
        return( c( PCA$values[1:k], PCA$vectors[,1:k] ) )
    }
    eigen.mat <- do.call( rbind, this.lapply( 1:floor(nrow(data)/win), .local ) )
    return( list(
            values = eigen.mat[1:k,],
            vectors = lapply(1:k,function(j) { eigen.mat[ncol(data)*k+(1:ncol(data)),] })
        ) )
}
