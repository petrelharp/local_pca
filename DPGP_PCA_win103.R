#!/home/rcf-40/hli465/bin/Rscript
#PBS -q cmb
#PBS -l walltime=100:00:00
#PBS -e /home/cmb-11/plr/hli465/Dpgp/all_samples_Chr3L/output
#PBS -o /home/cmb-11/plr/hli465/Dpgp/all_samples_Chr3L/output
#PBS -l nodes=1:ppn=1
#PBS -l mem=100gb,pmem=100gb,vmem=100gb
win <- 10^3
coded <- read.table("/home/cmb-11/plr/hli465/Dpgp/all_samples_Chr3L/coded_data_for_all_samples_seqs_both_low_NAs_Chr3L_with_SNP_Pos.txt",header=TRUE)
k <- 1:(floor(nrow(coded)/win))
usedata <- coded
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
PCs <- sapply(k, get.eigenvector, d=usedata)
write.table(PCs,"/home/cmb-11/plr/hli465/Dpgp/all_samples_Chr3L/fluctuation_PCA_win_103_Chr3L.txt",sep="\t")
