#!/home/rcf-40/hli465/bin/Rscript
#PBS -q cmb
#PBS -l walltime=200:00:00
#PBS -e /home/cmb-11/plr/hli465/FormatVCFbgz/output
#PBS -o /home/cmb-11/plr/hli465/FormatVCFbgz/output
#PBS -l nodes=1:ppn=1
#PBS -l mem=100gb,pmem=100gb,vmem=100gb
win <- 10^4
chr=read.table("/home/cmb-11/plr/hli465/FormatVCFbgz/chr1")
count=floor(nrow(chr)/win)
usedata <- chr[,-(1:2)]
get.eigenvector <- function(x, d) {
    chunk <- d[((x-1)*win + 1):(x*win), ]
    temp<-chunk
    temp<-data.matrix(temp)
    n<-ncol(temp)/2
    coded<-temp[,2*(1:n)]+temp[,2*(1:n)-1]
    data=coded
    M=rowMeans(data,na.rm=TRUE)
    M=rep(M,times=ncol(data))
    M=matrix(M,nrow=nrow(data),ncol=ncol(data),byrow=FALSE)
    data=data-M
    cov=cov(data,use="pairwise")
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
PCs <- sapply(1:count, get.eigenvector, d=usedata)
write.table(PCs,"/home/cmb-11/plr/hli465/FormatVCFbgz/fluctuation_PCA_win_104_chr1.txt",sep="\t")
