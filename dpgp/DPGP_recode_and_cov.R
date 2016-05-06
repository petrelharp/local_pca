#!/home/rcf-40/hli465/bin/Rscript
#PBS -q cmb
#PBS -l walltime=100:00:00
#PBS -e /home/cmb-11/plr/hli465/Dpgp/all_samples_Chr3L/output
#PBS -o /home/cmb-11/plr/hli465/Dpgp/all_samples_Chr3L/output
#PBS -l nodes=1:ppn=1
#PBS -l mem=100gb,pmem=100gb,vmem=100gb
data=read.table("/home/cmb-11/plr/hli465/Dpgp/all_samples_Chr3L/all_sample_seqs_Chr3L_with_SNP_Pos",header=TRUE,row.names=1)
geno <- as.matrix(data)
geno[ ! (geno %in% c("A","C","G","T")) ] <- NA
Gs <- rowSums(geno=="G",na.rm=TRUE)
As <- rowSums(geno=="A",na.rm=TRUE)
Cs <- rowSums(geno=="C",na.rm=TRUE)
Ts <- rowSums(geno=="T",na.rm=TRUE)
triallelic <- ( ( (Gs==0) + (As==0) + (Cs==0) + (Ts==0) ) > 2 )
maxcounts <- pmax(Gs,As,Cs,Ts)
major <- rep(NA,nrow(geno))
major[ Gs==maxcounts ] <- "G"
major[ As==maxcounts ] <- "A"
major[ Ts==maxcounts ] <- "T"
major[ Cs==maxcounts ] <- "C"
stopifnot(sum(is.na(major))==0)
coded <- matrix( NA, nrow=nrow(geno),ncol= ncol(geno) )
coded[ geno==major ] <- 0
coded[ geno!=major ] <- 1
coded[is.na(geno)]<-NA
colnames(coded)<-colnames(data)
rownames(coded)<-rownames(data)
write.table(coded,"/home/cmb-11/plr/hli465/Dpgp/all_samples_Chr3L/coded_data_all_samples_Chr3L_with_SNP_Pos.txt",sep="\t")
coded=data.matrix(coded)
b=apply(coded,2,function(x) sum(is.na(x))/length(x))
# remove individuals with more than this much missing data
kk=(b<0.08)
coded=coded[,kk]
a=apply(coded,1,function(x) sum(is.na(x))/length(x))
# remove sites with more than this much missing data
k=(a<0.2)
coded=coded[k,]
write.table(coded,"/home/cmb-11/plr/hli465/Dpgp/all_samples_Chr3L/coded_data_for_all_samples_seqs_both_low_NAs_Chr3L_with_SNP_Pos.txt",sep="\t")
data2=as.matrix(coded)
M=rowMeans(data2,na.rm=TRUE)
M=rep(M,times=ncol(data2))
M=matrix(M,nrow=nrow(data2),ncol=ncol(data2),byrow=FALSE)
data2=data2-M
cov=cov(data2,use="pairwise")
colnames(cov)<-colnames(coded)
write.table(cov,"/home/cmb-11/plr/hli465/Dpgp/all_samples_Chr3L/cov_data_for_all_samples_seqs_both_low_NAs_Chr3L.txt",sep="\t")

