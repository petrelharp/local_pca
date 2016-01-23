#!/home/rcf-40/hli465/bin/Rscript
#PBS -q cmb
#PBS -l walltime=20:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=100gb,pmem=100gb,vmem=100gb
chr<-read.table("/home/cmb-11/plr/hli465/FormatVCFbgz/chr1",na.strings="NA",stringsAsFactors=FALSE)
temp<-chr[,-(1:2)]
temp=data.matrix(temp)
n<-ncol(temp)/2
coded<-temp[,2*(1:n)]+temp[,2*(1:n)-1]
data=coded
M=rowMeans(data,na.rm=TRUE)
M=rep(M,times=ncol(data))
M=matrix(M,nrow=nrow(data),ncol=ncol(data),byrow=FALSE)
data=data-M
cov=cov(data,use="pairwise")
write.table(cov,"/home/cmb-11/plr/hli465/FormatVCFbgz/cov_data_for_chr1_new.txt",sep="\t")
