#!/home/rcf-40/hli465/bin/Rscript
#PBS -q cmb
#PBS -l walltime=100:00:00
#PBS -e /home/cmb-11/plr/hli465/POPRES/output
#PBS -o /home/cmb-11/plr/hli465/POPRES/output
#PBS -l nodes=1:ppn=1
#PBS -l mem=20gb,pmem=20gb,vmem=20gb
win <- 100
coded <- read.table("/home/cmb-11/plr/hli465/POPRES/coded_data_chr1.txt")
k <- 1:(floor(nrow(coded)/win))

get.eigenvector <- function(x, d) {
    step <- round(nrow(d)/10)
    chunk <- d[-(((x-1)*step + 1):(x*step)), ]
    temp<-chunk
    temp<-data.matrix(temp)
    data=temp
    M=rowMeans(data,na.rm=TRUE)
    M=rep(M,times=ncol(data))
    M=matrix(M,nrow=nrow(data),ncol=ncol(data),byrow=FALSE)
    data=data-M
    cov=cov(data,use="pairwise")
    if(sum(is.na(cov))>0) {return(rep(NA,nrow(cov)))}
    PCA=eigen(cov)
    Vec=PCA$vectors
    lam=PCA$values
    PC1=Vec[,1]
    return(PC1)
}

varfunction<-function(data){
    var=9/10*sum((data-mean(data))^2)
    return(var)
}

get.PC1s <- function(i, da) {

    index <- (100 * (i-1) + 1):(100 * i)
    usedata <- da[index, ]
    PC1s <- sapply(1:10, get.eigenvector, d=usedata)
    if(sum(is.na(PC1s))>0) {return(NA)}
    a=as.matrix(PC1s)
    a=t(a)
   A=matrix(0,nrow=nrow(a),ncol=ncol(a))
   S=rep(0,nrow(A))
  for(i in 1:nrow(A))
  {
        if (sum((a[1,]-a[i,])^2)<sum((a[1,]+a[i,])^2))
        {
                S[i]=1
                A[i,]=a[i,]     
        }
        else 
        {
                S[i]=-1
                A[i,]=-a[i,]    
        }
  }
     b=apply(A,2,varfunction)
    return(mean(b))
}

mbs <- sapply(k, get.PC1s, da=coded)
mbs=as.matrix(mbs)
rownames(mbs)=k
write.table(mbs,"/home/cmb-11/plr/hli465/POPRES/chr1_PC1_SE_win100_all",sep="\t")
