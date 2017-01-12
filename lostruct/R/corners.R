#' Find Points in a 2D Cloud Closest to the "Corners"
#'
#' Finds the k "extreme points" that lie closest to the minimal enclosing circle,
#' and returns a matrix with k columns containing the portion of points 
#' that lie closest to those extreme points.
#' If k is larger than the number of points on the minimal enclosing circle,
#' adds the closest points to that circle that are not in previously added points' neighborhoods
#' until k is reached.
#'
#' @param xy A two-column numeric matrix of coordinates.
#' @param prop The proportion of points to return for each corner.
#' @return A three column integer matrix giving the indices of the corresponding points.
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
