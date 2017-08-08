setwd("/home/cmb-12/fs3/hli465/Dpgp/simulations/")
win <- 10^4
library(lostruct)
coded <- vcf_windower("/home/rcf-40/pralph/panfs/local_pca/sims/neutral/hot_recomb.bcf",size=10^4,type='snp')
eigenstuff2 <- eigen_windows(coded, win=10^4, k=2,recode=FALSE)
windist2 <- pc_dist(eigenstuff2, npc=2)
write.table(windist2,file="/home/cmb-12/fs3/hli465/Dpgp/simulations/pairwise_distance_win104_hot_PCAk2",sep="\t")
eigenstuff5 <- eigen_windows(coded, win=10^4, k=5,recode=FALSE)
windist5 <- pc_dist(eigenstuff5, npc=5)
write.table(windist5,file="/home/cmb-12/fs3/hli465/Dpgp/simulations/pairwise_distance_win104_hot_PCAk5",sep="\t")

a<-as.matrix(windist2)
diag(a)=rep(0,nrow(a))                                           
a=sqrt(a)
windows_hotk2=rownames(a) 
fit2<-cmdscale(a,eig=TRUE, k=2)  
x2 <- fit2$points[,1]
y2 <- fit2$points[,2]
coordinate1_hotk2<- x2
coordinate1_hotk2=as.matrix(coordinate1_hotk2)
row.names(coordinate1_hotk2)=rownames(a)
coordinate2_hotk2<- y2
coordinate2_hotk2=as.matrix(coordinate2_hotk2)
row.names(coordinate2_hotk2)=rownames(a)

a<-as.matrix(windist5)
diag(a)=rep(0,nrow(a))                                           
a=sqrt(a)
windows_hotk5=rownames(a) 
fit2<-cmdscale(a,eig=TRUE, k=2)  
x2 <- fit2$points[,1]
y2 <- fit2$points[,2]
coordinate1_hotk5<- x2
coordinate1_hotk5=as.matrix(coordinate1_hotk5)
row.names(coordinate1_hotk2)=rownames(a)
coordinate2_hotk5<- y2
coordinate2_hotk5=as.matrix(coordinate2_hotk5)
row.names(coordinate2_hotk5)=rownames(a)


read.table("snp_pos_hot",header=TRUE)[,]->snp_pos
floor(length(snp_pos)/10^4) -> n
sapply(1:n, function(i){(snp_pos[i*10^4]+snp_pos[(i-1)*10^4+1])/2}) -> mid_pos
floor(mid_pos) -> wid_mid_pos #all wins
rwid_mid_pos_hotk2=as.matrix(wid_mid_pos[as.numeric(windows_hotk2)])/10^6     #wins really exsiting in "a"
rwid_mid_pos_hotk5=as.matrix(wid_mid_pos[as.numeric(windows_hotk5)])/10^6

pdf(file="MDS_plot_hot_sim_pos.pdf",width=10,height=6)
this.cex <- 0.6
layout(matrix(c(1:6), nrow=2,ncol=3,byrow=TRUE))
par(mar=c(4,4.5,2.5,2))
plot(coordinate1_hotk2, coordinate2_hotk2, xlab="coordinate1", ylab="coordinate2",pch=20,cex=this.cex)
plot(rwid_mid_pos_hotk2,coordinate1_hotk2,xlab="position (Mbp)",ylab="coordinate1",main="hot simulation PCAk2",pch=20,cex=this.cex)
plot(rwid_mid_pos_hotk2,coordinate2_hotk2,xlab="position (Mbp)",ylab="coordinate2",pch=20,cex=this.cex)

plot(coordinate1_hotk5, coordinate2_hotk5, xlab="coordinate1", ylab="coordinate2",pch=20,cex=this.cex)
plot(rwid_mid_pos_hotk5,coordinate1_hotk5,xlab="position (Mbp)",ylab="coordinate1",main="hot simulation PCAk5",pch=20,cex=this.cex)
plot(rwid_mid_pos_hotk5,coordinate2_hotk5,xlab="position (Mbp)",ylab="coordinate2",pch=20,cex=this.cex)

dev.off()


#delete windows within or around centromere (windows 426—522)
a<-as.matrix(windist2)
diag(a)=rep(0,nrow(a))                                           
a=sqrt(a)
a=a[-(426:522),-(426:522)]
windows_hotk2=rownames(a) 
fit2<-cmdscale(a,eig=TRUE, k=2)  
x2 <- fit2$points[,1]
y2 <- fit2$points[,2]
coordinate1_hotk2<- x2
coordinate1_hotk2=as.matrix(coordinate1_hotk2)
row.names(coordinate1_hotk2)=rownames(a)
coordinate2_hotk2<- y2
coordinate2_hotk2=as.matrix(coordinate2_hotk2)
row.names(coordinate2_hotk2)=rownames(a)

a<-as.matrix(windist5)
diag(a)=rep(0,nrow(a))                                           
a=sqrt(a)
a=a[-(426:522),-(426:522)]
windows_hotk5=rownames(a) 
fit2<-cmdscale(a,eig=TRUE, k=2)  
x2 <- fit2$points[,1]
y2 <- fit2$points[,2]
coordinate1_hotk5<- x2
coordinate1_hotk5=as.matrix(coordinate1_hotk5)
row.names(coordinate1_hotk2)=rownames(a)
coordinate2_hotk5<- y2
coordinate2_hotk5=as.matrix(coordinate2_hotk5)
row.names(coordinate2_hotk5)=rownames(a)


read.table("~/Documents/Drosophila/sims/snp_pos_hot",header=TRUE)[,]->snp_pos
floor(length(snp_pos)/10^4) -> n
sapply(1:n, function(i){(snp_pos[i*10^4]+snp_pos[(i-1)*10^4+1])/2}) -> mid_pos
floor(mid_pos) -> wid_mid_pos #all wins
rwid_mid_pos_hotk2=as.matrix(wid_mid_pos[as.numeric(windows_hotk2)])/10^6     #wins really exsiting in "a"
rwid_mid_pos_hotk5=as.matrix(wid_mid_pos[as.numeric(windows_hotk5)])/10^6

pdf(file="MDS_plot_hot_sim_without_centromere_pos.pdf",width=10,height=6)
this.cex <- 0.6
layout(matrix(c(1:6), nrow=2,ncol=3,byrow=TRUE))
par(mar=c(4,4.5,2.5,2))
plot(coordinate1_hotk2, coordinate2_hotk2, xlab="coordinate1", ylab="coordinate2",pch=20,cex=this.cex)
plot(rwid_mid_pos_hotk2,coordinate1_hotk2,xlab="position (Mbp)",ylab="coordinate1",main="hot simulation PCAk2",pch=20,cex=this.cex)
plot(rwid_mid_pos_hotk2,coordinate2_hotk2,xlab="position (Mbp)",ylab="coordinate2",pch=20,cex=this.cex)

plot(coordinate1_hotk5, coordinate2_hotk5, xlab="coordinate1", ylab="coordinate2",pch=20,cex=this.cex)
plot(rwid_mid_pos_hotk5,coordinate1_hotk5,xlab="position (Mbp)",ylab="coordinate1",main="hot simulation PCAk5",pch=20,cex=this.cex)
plot(rwid_mid_pos_hotk5,coordinate2_hotk5,xlab="position (Mbp)",ylab="coordinate2",pch=20,cex=this.cex)

dev.off()



