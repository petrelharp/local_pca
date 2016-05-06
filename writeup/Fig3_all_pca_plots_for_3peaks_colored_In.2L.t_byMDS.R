setwd("~/Documents/Drosophila/Chr2L")
c1=as.matrix(read.table("Chr2L_recomchunk_win103_cov1_byMDS_ordered.txt"))
c2=as.matrix(read.table("Chr2L_recomchunk_win103_cov2_byMDS_ordered.txt"))
c3=as.matrix(read.table("Chr2L_recomchunk_win103_cov3_byMDS_ordered.txt"))
#c=as.matrix(read.table("cov_data_for_all_samples_seqs_both_low_NAs_Chr2L.txt"))

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
#PCA=eigen(c)
#Vec=PCA$vectors
#lam=PCA$values
#PC1=as.matrix(Vec[,1])
#PC2=as.matrix(Vec[,2])

origin=as.matrix(colnames(c1))
dir=as.matrix(read.table("~/Documents/Drosophila/samples_inversions.txt"))
origin1=unlist(strsplit(origin[,1], split="_"))
origin2=origin1[2*(1:nrow(c1))-1]
inv1=match(origin2,dir[,1])
inv2=dir[inv1,2]
library(RColorBrewer)
display.brewer.all(3,colorblindFriendly=TRUE)
cols <- brewer.pal(6,"Set2")[4:6]
pdf(file="Fig3_all_pca_plots_for_Chr2L_3peaks_color_by_In.2L.t_byMDS.pdf",width=12,height=4)
layout(matrix(c(1,2,3), nrow=1,byrow=TRUE))
group=as.numeric(as.factor(inv2))
group[which(is.na(group))]=3
par(mar=c(5,4,3,2))
plot(PC11,PC12,pch=group,col= cols[group])
plot(PC21,PC22,pch=group,col= cols[group],main="Drosophila chromosome 2L")
plot(PC31,PC32,pch=group,col= cols[group])
legend(0.048,0.09,pch=1:3,col= cols,legend=c("INV","ST","N"))
dev.off()


