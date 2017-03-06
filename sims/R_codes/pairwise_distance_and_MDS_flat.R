setwd("/home/cmb-12/fs3/hli465/Dpgp/simulations/")
win <- 10^4
library(lostruct)
coded <- vcf_windower("/home/rcf-40/pralph/panfs/local_pca/sims/neutral/flat_recomb.bcf",size=10^4,type='snp')
eigenstuff2 <- eigen_windows(coded, win=10^4, k=2,recode=FALSE)
windist2 <- pc_dist(eigenstuff2, npc=2)
write.table(windist2,file="/home/cmb-12/fs3/hli465/Dpgp/simulations/pairwise_distance_win104_flat_PCAk2",sep="\t")
eigenstuff5 <- eigen_windows(coded, win=10^4, k=5,recode=FALSE)
windist5 <- pc_dist(eigenstuff5, npc=5)
write.table(windist5,file="/home/cmb-12/fs3/hli465/Dpgp/simulations/pairwise_distance_win104_flat_PCAk5",sep="\t")

a<-as.matrix(windist2)
diag(a)=rep(0,nrow(a))                                           
a=sqrt(a)
windows_flatk2=rownames(a) 
fit2<-cmdscale(a,eig=TRUE, k=2)  
x2 <- fit2$points[,1]
y2 <- fit2$points[,2]
coordinate1_flatk2<- x2
coordinate1_flatk2=as.matrix(coordinate1_flatk2)
row.names(coordinate1_flatk2)=rownames(a)
coordinate2_flatk2<- y2
coordinate2_flatk2=as.matrix(coordinate2_flatk2)
row.names(coordinate2_flatk2)=rownames(a)

a<-as.matrix(windist5)
diag(a)=rep(0,nrow(a))                                           
a=sqrt(a)
windows_flatk5=rownames(a) 
fit2<-cmdscale(a,eig=TRUE, k=2)  
x2 <- fit2$points[,1]
y2 <- fit2$points[,2]
coordinate1_flatk5<- x2
coordinate1_flatk5=as.matrix(coordinate1_flatk5)
row.names(coordinate1_flatk2)=rownames(a)
coordinate2_flatk5<- y2
coordinate2_flatk5=as.matrix(coordinate2_flatk5)
row.names(coordinate2_flatk5)=rownames(a)


read.table("snp_pos_flat",header=TRUE)[,]->snp_pos
floor(length(snp_pos)/10^4) -> n
sapply(1:n, function(i){(snp_pos[i*10^4]+snp_pos[(i-1)*10^4+1])/2}) -> mid_pos
floor(mid_pos) -> wid_mid_pos #all wins
rwid_mid_pos_flatk2=as.matrix(wid_mid_pos[as.numeric(windows_flatk2)])/10^6     #wins really exsiting in "a"
rwid_mid_pos_flatk5=as.matrix(wid_mid_pos[as.numeric(windows_flatk5)])/10^6

pdf(file="MDS_plot_flat_sim_pos.pdf",width=10,height=6)
this.cex <- 0.6
layout(matrix(c(1:6), nrow=2,ncol=3,byrow=TRUE))
par(mar=c(4,4.5,2.5,2))
plot(coordinate1_flatk2, coordinate2_flatk2, xlab="coordinate1", ylab="coordinate2",pch=20,cex=this.cex)
plot(rwid_mid_pos_flatk2,coordinate1_flatk2,xlab="position (Mbp)",ylab="coordinate1",main="flat simulation PCAk2",pch=20,cex=this.cex)
plot(rwid_mid_pos_flatk2,coordinate2_flatk2,xlab="position (Mbp)",ylab="coordinate2",pch=20,cex=this.cex)

plot(coordinate1_flatk5, coordinate2_flatk5, xlab="coordinate1", ylab="coordinate2",pch=20,cex=this.cex)
plot(rwid_mid_pos_flatk5,coordinate1_flatk5,xlab="position (Mbp)",ylab="coordinate1",main="flat simulation PCAk5",pch=20,cex=this.cex)
plot(rwid_mid_pos_flatk5,coordinate2_flatk5,xlab="position (Mbp)",ylab="coordinate2",pch=20,cex=this.cex)

dev.off()
