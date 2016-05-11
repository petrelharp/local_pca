#' Find Matrix of Distances Between Matrices from Eigenvector/values
#'
#' Given a matrix of sets of eigenvalues/vectors for a set of n symmetric matrices as output by eigen_windows(),
#' return the (n x n) matrix whose [i,j]th element is the Frobenius norm of the difference between the i-th and the j-th matrices,
#' as approximated by the given eigenvalues/vectors. The approximation is in pseudocode
#'      M = values[1] * outer(vectors[1],vectors[1]) + values[2] * outer(vectors[2],vectors[2]) + ...
#' If \code{normalize} is TRUE, then the the matrices are normalized to have norm 1, so that in the definition of M,
#' \code{values} is replaced by \code{values/sqrt(sum(values^2))}.
#'
#' @param x Matrix as output by eigen_windows().
#' @param npc Number of pcas computed.
#' @param normalize Normalize the matrices to have the same norm?  Can be "L1" (default), "L2", or FALSE.
#' @param mc.cores If this is greater than 1, parallel::mclapply will be used.
#' @return A symmetric, numeric matrix with number of columns equal to the number of columns in eigen.win$values.
#' @export
pc_dist <- function( x, npc, normalize="L1", mc.cores=1 ) {
    this.lapply <- if (mc.cores>1) { function (...) parallel::mclapply(...,mc.cores=mc.cores) } else { lapply }
    values <- x[,1:npc,drop=FALSE]
    vectors <- x[,-(1:npc)]
    if (normalize=="L1") { values <- sweep( values, 1, abs(rowSums(values)), "/" ) }
    if (normalize=="L2") { values <- sweep( values, 1, sqrt(rowSums(values^2)), "/" ) }
    n <- nrow(values)
    # could speed this up by a factor of two taking advantage of symmetry
    emat <- function (u) { matrix(u,ncol=npc) }
    out <- do.call( rbind, this.lapply( 1:n, function (i) {
                sapply( 1:n, function (j) {
                    dist_from_pcs( values[i,], emat(vectors[i,]), values[j,], emat(vectors[j,]) )
                } )
            } ) )
    # symmetrize
    return( (out + t(out))/2 )
}


#' Suppose that (u) and (v) are sets of vectors, and that 
#'    A = a_1 * u_1 u_1^T + ... + a_j * u_j u_j^T = U diag(a) U^T
#'    B = b_1 * v_1 v_1^T + ... + b_j * v_k v_k^T = V diag(b) V^T,
#' where U is the matrix whose columns are (u), and likewise for V.
#' Then 
#'    (A-B)^T (A-B) = A^T A + B^T B - A^T B - B^T A
#' and so
#'    ||A-B||^2 = ||A||^2 + ||B||^2 - 2 tr( A^T B )
#' By the cyclic invariance of trace, if X = U^T V then
#'    tr( A^T B ) = tr( U diag(a) U^T V diag(b) V^T ) = tr( diag(a) X diag(b) X^T )
#'
#' The code assumes that vectors1 and vectors2 are each orthonormal.
dist_from_pcs <- function (values1,vectors1,values2,vectors2) {
    # this is diag(sqrt(b)) X^T
    bXt <- sqrt(values2) * crossprod(vectors2,vectors1)
    return( sum(values1^2) + sum(values2^2) - 2 * sum( values1 * colSums(bXt^2) ) )
}
