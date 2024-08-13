#' Find Points in a 2D Cloud Closest to the "Corners"
#'
#' Finds the 'k' "extreme points" (the "corners") that lie closest
#' to the minimal circle enclosing the input collection of points
#' in two dimensions, and returns for each of these corners
#' which of the input points are closest to that corner.
#' The return value is a matrix with 'k' columns,
#' whose 'i'-th column contains the indices of those points
#' that lie closest to the 'i'-th corner.
#'
#' If the output is 'out', then 'out[,i]' is a vector of
#' length 'nrow(xy) * prop' for which 'xy[ out[,i], ]' is
#' the coordinates of the closest points to the 'i'th corner.
#' The rows of 'out' are in no particular order.
#'
#' If k is larger than the number of points on the minimal enclosing circle,
#' this adds the closest points to that circle that are not in previously
#' added points' neighborhoods until k is reached.
#' If k is less than 3, this will behave as if k=3.
#'
#' @param xy A two-column numeric matrix of coordinates.
#' @param prop The proportion of points to return for each corner.
#' @param k The number of corners (should be at least 3).
#' @return A k-column integer matrix whose i-th column gives the indices of the points closest to the i-th corner.
#'
#' @examples
#' xy <- runif(200) |> matrix(ncol=2)
#' out <- corners(xy, k=4, prop=0.05)
#' \dontrun{
#' plot(xy)
#' for (i in 1:4) points(xy[out[,i],], col=i, pch=20, cex=2)
#' }
#' @export
corners <- function (xy, prop, k=3) {
    xy <- as.matrix(xy)
    mincirc <- enclosing_circle( xy )
    cidx <- mincirc$index
    dpt <- function (uv) { sqrt( (xy[,1]-uv[1])^2 + (xy[,2]-uv[2])^2 ) }
    closest <- function (u) { which( (!is.na(u)) & (u <= quantile(u,probs=prop,na.rm=TRUE)) ) }
    dists <- sapply( cidx, function (k) { dpt(xy[k,]) } )
    out <- apply( dists, 2, closest )
    if (length(cidx)<k) {
        cdists <- dpt( mincirc$ctr )
        use.these <- setdiff( 1:nrow(xy), unlist(out) )
    }
    while ( length(cidx)<k && any(use.these) ) {
        new.idx <- use.these[ which.max( cdists[use.these] ) ]
        cidx <- c( cidx, new.idx )
        new.dists <- dpt( xy[new.idx,] )
        out <- cbind( out, closest(new.dists) )
        use.these <- setdiff( use.these, out[,ncol(out)] )
    }
    return(out)
}
