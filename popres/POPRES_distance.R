#!/home/rcf-40/hli465/bin/Rscript
#PBS -q cmb
#PBS -l walltime=100:00:00
#PBS -e /home/cmb-11/plr/hli465/POPRES/fluctuation_PCA/output
#PBS -o /home/cmb-11/plr/hli465/POPRES/fluctuation_PCA/output
#PBS -l nodes=1:ppn=1
#PBS -l mem=100gb,pmem=100gb,vmem=100gb
PCAs=as.matrix(read.table("fluctuation_PCA_win_100_chr1.txt"))
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
write.table(Distance,file="quick_method_pairwise_distance_between_win_100_chr1",sep="\t")
