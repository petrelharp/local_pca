#!/home/rcf-40/hli465/bin/Rscript
#PBS -q cmb
#PBS -l walltime=600:00:00
#PBS -e /home/cmb-11/plr/hli465/FormatVCFbgz/fluctuation_PCA/output
#PBS -o /home/cmb-11/plr/hli465/FormatVCFbgz/fluctuation_PCA/output
#PBS -l nodes=1:ppn=1
#PBS -l mem=100gb,pmem=100gb,vmem=100gb
PCA_1=as.matrix(read.table("/home/cmb-11/plr/hli465/FormatVCFbgz/fluctuation_PCA/fluctuation_PCA_win_104_chr1.txt"))
PCA_2=as.matrix(read.table("/home/cmb-11/plr/hli465/FormatVCFbgz/fluctuation_PCA/fluctuation_PCA_win_104_chr2.txt"))
PCA_3=as.matrix(read.table("/home/cmb-11/plr/hli465/FormatVCFbgz/fluctuation_PCA/fluctuation_PCA_win_104_chr3.txt"))
PCA_4=as.matrix(read.table("/home/cmb-11/plr/hli465/FormatVCFbgz/fluctuation_PCA/fluctuation_PCA_win_104_chr4.txt"))
PCA_5=as.matrix(read.table("/home/cmb-11/plr/hli465/FormatVCFbgz/fluctuation_PCA/fluctuation_PCA_win_104_chr5.txt"))
PCA_6=as.matrix(read.table("/home/cmb-11/plr/hli465/FormatVCFbgz/fluctuation_PCA/fluctuation_PCA_win_104_chr6.txt"))
PCA_7=as.matrix(read.table("/home/cmb-11/plr/hli465/FormatVCFbgz/fluctuation_PCA/fluctuation_PCA_win_104_chr7.txt"))
PCA_8=as.matrix(read.table("/home/cmb-11/plr/hli465/FormatVCFbgz/fluctuation_PCA/fluctuation_PCA_win_104_chr8.txt"))
PCAs=cbind(PCA_1, PCA_2, PCA_3, PCA_4, PCA_5, PCA_6, PCA_7, PCA_8)
index <- (1:ncol(PCAs))[!is.na(PCAs[1,])]
k=which(is.na(PCAs[1,]))
if(length(k)==0) {PCA=PCAs}
temp=nrow(PCA)/2
PC1=PCA[1:(temp-1),]
PC2=PCA[(temp+1):(2*temp-1),]
Lam=rbind(PCA[temp,],PCA[2*temp,])
count=ncol(PC1)
Distance=matrix(NA,nrow=count,ncol=count)
for(i in 1:count)
{
        V1=sqrt(Lam[1,i]/(Lam[1,i]+Lam[2,i]))*PC1[,i]
        V2=sqrt(Lam[2,i]/(Lam[1,i]+Lam[2,i]))*PC2[,i]
        part1=(sum(V1*V1))^2+(sum(V2*V2))^2
        for(j in i:count)
        {
                U1=sqrt(Lam[1,j]/(Lam[1,j]+Lam[2,j]))*PC1[,j]
                U2=sqrt(Lam[2,j]/(Lam[1,j]+Lam[2,j]))*PC2[,j]
                part2=(sum(U1*U1))^2+(sum(U2*U2))^2
                part3=(sum(V1*U1))^2+(sum(V1*U2))^2+(sum(V2*U1))^2+(sum(V2*U2))^2
                Distance[i,j]=part1+part2-2*part3
        }
}
write.table(Distance,file="/home/cmb-11/plr/hli465/FormatVCFbgz/fluctuation_PCA/quick_method_pairwise_distance_between_win_104_all_chr_medicago",sep="\t")
