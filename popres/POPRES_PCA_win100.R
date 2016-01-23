#!/home/rcf-40/hli465/bin/Rscript
#PBS -q cmb
#PBS -l walltime=100:00:00
#PBS -e /home/cmb-11/plr/hli465/POPRES/output
#PBS -o /home/cmb-11/plr/hli465/POPRES/output
#PBS -l nodes=1:ppn=1
#PBS -l mem=20gb,pmem=20gb,vmem=20gb

win <- 100
coded <- read.table("coded_data_chr1.txt")
k <- 1:(floor(nrow(coded)/win))

#' Given the index of a window and the data frame,
#' returns the first two eigenvectors and eigenvalues of the covariance matrix of that window
#' in the order (vec1,value1,vec2,value2).
#'
get.eigenvector <- function(x, d) {
    chunk <- d[((x-1)*win + 1):(x*win), ]
    temp<-chunk
    temp<-data.matrix(temp)
    data=temp
    M=rowMeans(data,na.rm=TRUE)
    M=rep(M,times=ncol(data))
    M=matrix(M,nrow=nrow(data),ncol=ncol(data),byrow=FALSE)
    data=data-M
    cov=cov(data,use="pairwise")
    if(sum(is.na(cov))>0) {return(rep(NA,2*(nrow(cov)+1)))}
    PCA=eigen(cov)
    Vec=PCA$vectors
    lam=PCA$values
    PC1=Vec[,1]
    PC2=Vec[,2]
    lam1=lam[1]
    lam2=lam[2]
    PCs=c(PC1,lam1,PC2,lam2)
    return(PCs)
}

PCs <- sapply(k, get.eigenvector, d=coded)
write.table(PCs,"fluctuation_PCA_win_100_chr1.txt",sep="\t")
