#' Find Points in a 2D Cloud Closest to the "Corners"
#'
#' Finds the three "extreme points" that lie on the minimal enclosing circle,
#' and returns a list of length three containing the portion of points 
#' that lie closest to those extreme points.
#'
#' @param xy A two-column numeric matrix of coordinates.
#' @param prop The proportion of points to return for each corner.
#' @return A three column integer matrix giving the indices of the corresponding points.
#' @export
corners <- function (xy, prop) {
    mincirc <- enclosing_circle( xy )
    cidx <- mincirc$index
    dists <- sqrt( outer( xy[,1], xy[cidx,1], "-" )^2 + outer( xy[,2], xy[cidx,2], "-" )^2 )
    return( apply( dists, 2, function (u) {
                which( u <= quantile(u,probs=prop) )
            } ) )
}

