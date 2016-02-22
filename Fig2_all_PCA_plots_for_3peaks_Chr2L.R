setwd("~/Documents/Drosophila/Chr2L")
c1=as.matrix(read.table("Chr2L_recomchunk_win103_cov1.txt"))
c2=as.matrix(read.table("Chr2L_recomchunk_win103_cov2.txt"))
c3=as.matrix(read.table("Chr2L_recomchunk_win103_cov3.txt"))
c=as.matrix(read.table("cov_data_for_all_samples_seqs_both_low_NAs_Chr2L.txt"))

PCA1=eigen(c1)
Vec1=PCA1$vectors
lam1=PCA1$values
PC11=as.matrix(Vec1[,1])
PC12=as.matrix(Vec1[,2])
PCA2=eigen(c2)
Vec2=PCA2$vectors
lam2=PCA2$values
PC21=as.matrix(Vec2[,1])
PC22=as.matrix(Vec2[,2])
PCA3=eigen(c3)
Vec3=PCA3$vectors
lam3=PCA3$values
PC31=as.matrix(Vec3[,1])
PC32=as.matrix(Vec3[,2])

PCA=eigen(c)
Vec=PCA$vectors
lam=PCA$values
PC1=as.matrix(Vec[,1])
PC2=as.matrix(Vec[,2])

origin=colnames(c) 
origin1=substring(origin,1,2)
countrys=as.matrix(read.table("population_country.txt"))
for(i in 1:nrow(countrys)){origin1[which(origin1==countrys[i,1])]=countrys[i,2]}
rownames(PC1)=origin1
rownames(PC2)=origin1
pdf(file="Fig2_all_pca_plots_for_Chr2L_3peaks.pdf",width=20,height=16)
layout(matrix(c(1,2,3,4), nrow=2,byrow=TRUE))
group=as.numeric(as.factor(origin1))
par(mar=c(6,6,5,5),xpd=TRUE)
plot(PC11,PC12,pch=group,col=rainbow(16)[group])
plot(PC21,PC22,pch=group,col=rainbow(16)[group])
plot(PC31,PC32,pch=group,col=rainbow(16)[group])
plot(PC1,PC2,pch=group,col=rainbow(16)[group])
legend(0.08,0.15,pch=1:16,col=rainbow(16),legend=levels(factor(origin1)))
dev.off()













