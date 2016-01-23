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
#' @param mc.cores If this is greater than 1, parallel::mclapply will be used.
#' @return A named list, with components
#'   $values : A numeric matrix of eigenvalues, one column for each window and k rows.
#'   $vectors : A numeric matrix, one column for each windows, having the eigenvectors in order (so it has k * ncol(matrix) rows).
#' @export
eigen_windows <- function ( data, win, k, mc.cores=1 ) {
    this.lapply <- if (mc.cores>1) { function (...) parallel::mclapply(...,mc.cores=mc.cores) } else { lapply }
    data <- sweep( data, 1, rowMeans(data,na.rm=TRUE) )
    .local <- function(x) {
        chunk <- data[((x-1)*win + 1):(x*win), ]
        cov <- cov(chunk,use="pairwise")
        if(any(is.na(cov))) {return(rep(NA,2*(nrow(cov)+1)))}
        PCA <- eigen(cov)
        # returns in order (values, vectors)
        return( c( PCA$values[1:k], PCA$vectors[,1:k] ) )
    }
    eigen.mat <- do.call( cbind, this.lapply( 1:floor(nrow(data)/win), .local ) )
    return( list(
            values = eigen.mat[1:k,],
            vectors = eigen.mat[-(1:k),]
        ) )
}
