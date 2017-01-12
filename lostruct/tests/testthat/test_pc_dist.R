
context("Testing pc_dist")

set.seed(42)
n <- 10

x <- cov( matrix(rnorm(2*n^2), ncol=n) )
y <- cov( matrix(rnorm(2*n^2), ncol=n) )
x.eig <- eigen(x)
y.eig <- eigen(y)

elist <- function (k) {
    # transform to the format pc_dist wants
    rbind( c(sum(x^2), x.eig$values[1:k], as.vector( x.eig$vectors[,1:k] ) ),
           c(sum(y^2), y.eig$values[1:k], as.vector( y.eig$vectors[,1:k] ) ) )
}

f <- function (k,eig) {
    out <- eig$values[1] * outer(eig$vectors[,1],eig$vectors[,1],"*")
    for (k in seq_len(k)[-1]) {
        out <- out + eig$values[k] * outer(eig$vectors[,k],eig$vectors[,k],"*")
    }
    return(out)
}
g1 <- function (k,eig) {
    eig$values <- eig$values/sum(abs(eig$values[1:k]))
    f(k,eig)
}
g2 <- function (k,eig) {
    eig$values <- eig$values/sqrt(sum(eig$values[1:k]^2))
    f(k,eig)
}

kk <- c(1,5,10)
for (k in kk) {
    # non-normalized
    ff <- sum( (f(k,x.eig) - f(k,y.eig))^2 )
    fm <- matrix(c(0,ff,ff,0),nrow=2)
    expect_equal( fm, pc_dist(elist(k),npc=k,normalize=FALSE) )
    # normalized, L1
    gg <- sum( (g1(k,x.eig) - g1(k,y.eig))^2 )
    gm <- matrix(c(0,gg,gg,0),nrow=2)
    expect_equal( gm, pc_dist(elist(k),npc=k,normalize="L1") )
    # normalized, L2
    gg <- sum( (g2(k,x.eig) - g2(k,y.eig))^2 )
    gm <- matrix(c(0,gg,gg,0),nrow=2)
    expect_equal( gm, pc_dist(elist(k),npc=k,normalize="L2") )
}

## weighted

w <- rexp( n )
x.eig <- eigen(sqrt(w) * x %*% diag(sqrt(w)))
x.eig$vectors <- diag(1/sqrt(w)) %*% x.eig$vectors
y.eig <- eigen(sqrt(w) * y %*% diag(sqrt(w)))
y.eig$vectors <- diag(1/sqrt(w)) %*% y.eig$vectors

# check vectors orthogonal in l2(w)
expect_equal( crossprod( diag(sqrt(w)) %*% x.eig$vectors ), diag(ncol(x.eig$vectors)) )
expect_equal( crossprod( diag(sqrt(w)) %*% y.eig$vectors ), diag(ncol(y.eig$vectors)) )
# check have an expression for x, y
expect_equal( x, x.eig$vectors %*% ( x.eig$values * t( x.eig$vectors ) ) )
expect_equal( y, y.eig$vectors %*% ( y.eig$values * t( y.eig$vectors ) ) )
# check partition of sum of squares
expect_equal( sum( w * (x^2 %*% diag(w)) ), sum( x.eig$values^2 ) )
expect_equal( sum( w * (y^2 %*% diag(w)) ), sum( y.eig$values^2 ) )

kk <- c(1,5,10)
for (k in kk) {
    # non-normalized
    ff <- sum(  w * sweep( (f(k,x.eig) - f(k,y.eig))^2, 2, w, "*" ) )
    fm <- matrix(c(0,ff,ff,0),nrow=2)
    expect_equal( fm, pc_dist(elist(k),npc=k,w=w,normalize=FALSE) )
    # normalized, L1
    gg <- sum( w * sweep( (g1(k,x.eig) - g1(k,y.eig))^2, 2, w, "*" ) )
    gm <- matrix(c(0,gg,gg,0),nrow=2)
    expect_equal( gm, pc_dist(elist(k),npc=k,w=w,normalize="L1") )
    # normalized, L2
    gg <- sum( w * sweep( (g2(k,x.eig) - g2(k,y.eig))^2, 2, w, "*" ) )
    gm <- matrix(c(0,gg,gg,0),nrow=2)
    expect_equal( gm, pc_dist(elist(k),npc=k,w=w,normalize="L2") )
}

