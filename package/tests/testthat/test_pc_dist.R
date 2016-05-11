
context("Testing pc_dist")

set.seed(42)
n <- 10

x <- cov( matrix(rnorm(2*n^2), ncol=n) )
y <- cov( matrix(rnorm(2*n^2), ncol=n) )
x.eig <- eigen(x)
y.eig <- eigen(y)

elist <- function (k) {
    # transform to the format pc_dist wants
    rbind( c(x.eig$values[1:k], as.vector( x.eig$vectors[,1:k] ) ),
           c(y.eig$values[1:k], as.vector( y.eig$vectors[,1:k] ) ) )
}

f <- function (k,eig) {
    out <- eig$values[1] * outer(eig$vectors[,1],eig$vectors[,1],"*")
    for (k in seq_len(k)[-1]) {
        out <- out + eig$values[k] * outer(eig$vectors[,k],eig$vectors[,k],"*")
    }
    return(out)
}
g <- function (k,eig) {
    eig$values <- eig$values/sqrt(sum(eig$values[1:k]^2))
    f(k,eig)
}

kk <- c(1,5,10)
for (k in kk) {
    # non-normalized
    ff <- sum( (f(k,x.eig) - f(k,y.eig))^2 )
    fm <- matrix(c(0,ff,ff,0),nrow=2)
    expect_equal( fm, pc_dist(elist(k),npc=k,normalize=FALSE) )
    # normalized
    gg <- sum( (g(k,x.eig) - g(k,y.eig))^2 )
    gm <- matrix(c(0,gg,gg,0),nrow=2)
    expect_equal( gm, pc_dist(elist(k),npc=k,normalize=TRUE) )
}
