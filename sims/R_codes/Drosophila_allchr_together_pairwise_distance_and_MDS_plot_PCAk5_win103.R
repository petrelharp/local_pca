win <- 10^3
library(lostruct)
fn=c("Chr2L","Chr2R","Chr3L","Chr3R","ChrX")
for(i in 1:length(fn))
{
   coded <- read.table(gettextf("/home/cmb-12/fs3/hli465/Dpgp/all_samples_%s/coded_data_for_all_samples_seqs_both_low_NAs_%s_with_SNP_Pos.txt",fn[i],fn[i]),header=TRUE)   
   eigenstuff <- eigen_windows(coded, win=10^3, k=5,recode=FALSE)
   windist <- pc_dist(eigenstuff, npc=5)
   write.table(windist,file=gettextf("/home/cmb-12/fs3/hli465/Dpgp/all_samples_%s/pairwise_distance_between_win103_%s_PCAk5",fn[i],fn[i]),sep="\t")
}

#Chr2L
data<-read.table("/home/cmb-12/fs3/hli465/Dpgp/all_samples_Chr2L/pairwise_distance_between_win103_Chr2L_PCAk5",check.names=FALSE)
ind1 <- apply(data, 1, function(x) all(is.na(x)))
a1=data[!ind1,]
ind2 <- apply(a1, 2, function(x) all(is.na(x)))
a2=a1[,!ind2]  
diag(a2)=rep(0,nrow(a2))                                           #delete NAs
a=sqrt(a2)

a <- as.matrix(a)
windows=rownames(a) 
fit2<-cmdscale(a,eig=TRUE, k=2)  
x2 <- fit2$points[,1]
y2 <- fit2$points[,2]
coordinate1<- x2
coordinate1=as.matrix(coordinate1)
row.names(coordinate1)=rownames(a)
coordinate2<- y2
coordinate2=as.matrix(coordinate2)
row.names(coordinate2)=rownames(a)
coordinate1_2L=coordinate1
coordinate2_2L=coordinate2

read.table("/home/cmb-12/fs3/hli465/Dpgp/all_samples_Chr2L/SNP_Pos_Chr2L.txt")[[1]]->snp_pos
floor(length(snp_pos)/1000) -> n
sapply(1:n, function(i){(snp_pos[i*1000]+snp_pos[(i-1)*1000+1])/2}) -> mid_pos
floor(mid_pos) -> wid_mid_pos #all wins
rwid_mid_pos_2L=as.matrix(wid_mid_pos[as.numeric(windows)])/10^6       #wins really exsiting in "a"
row.names(rwid_mid_pos_2L)=rownames(a)


#Chr2R
data<-read.table("/home/cmb-12/fs3/hli465/Dpgp/all_samples_Chr2R/pairwise_distance_between_win103_Chr2R_PCAk5",check.names=FALSE)
ind1 <- apply(data, 1, function(x) all(is.na(x)))
a1=data[!ind1,]
ind2 <- apply(a1, 2, function(x) all(is.na(x)))
a2=a1[,!ind2]  
diag(a2)=rep(0,nrow(a2))                                           #delete NAs
a=sqrt(a2)

a <- as.matrix(a)
windows=rownames(a) 
fit2<-cmdscale(a,eig=TRUE, k=2)  
x2 <- fit2$points[,1]
y2 <- fit2$points[,2]
coordinate1<- x2
coordinate1=as.matrix(coordinate1)
row.names(coordinate1)=rownames(a)
coordinate2<- y2
coordinate2=as.matrix(coordinate2)
row.names(coordinate2)=rownames(a)
coordinate1_2R=coordinate1
coordinate2_2R=coordinate2
read.table("/home/cmb-12/fs3/hli465/Dpgp/all_samples_Chr2R/SNP_Pos_Chr2R.txt")[[1]]->snp_pos
floor(length(snp_pos)/1000) -> n
sapply(1:n, function(i){(snp_pos[i*1000]+snp_pos[(i-1)*1000+1])/2}) -> mid_pos
floor(mid_pos) -> wid_mid_pos #all wins
rwid_mid_pos_2R=as.matrix(wid_mid_pos[as.numeric(windows)])/10^6       #wins really exsiting in "a"
row.names(rwid_mid_pos_2R)=rownames(a)

#Chr3L
data<-read.table("/home/cmb-12/fs3/hli465/Dpgp/all_samples_Chr3L/pairwise_distance_between_win103_Chr3L_PCAk5",check.names=FALSE)
ind1 <- apply(data, 1, function(x) all(is.na(x)))
a1=data[!ind1,]
ind2 <- apply(a1, 2, function(x) all(is.na(x)))
a2=a1[,!ind2] 
diag(a2)=rep(0,nrow(a2))                                            #delete NAs
a=sqrt(a2)

a <- as.matrix(a)
windows=rownames(a) 
fit2<-cmdscale(a,eig=TRUE, k=2)  
x2 <- fit2$points[,1]
y2 <- fit2$points[,2]
coordinate1<- x2
coordinate1=as.matrix(coordinate1)
row.names(coordinate1)=rownames(a)
coordinate2<- y2
coordinate2=as.matrix(coordinate2)
row.names(coordinate2)=rownames(a)
coordinate1_3L=coordinate1
coordinate2_3L=coordinate2
read.table("/home/cmb-12/fs3/hli465/Dpgp/all_samples_Chr3L/SNP_Pos_Chr3L.txt")[[1]]->snp_pos
floor(length(snp_pos)/1000) -> n
sapply(1:n, function(i){(snp_pos[i*1000]+snp_pos[(i-1)*1000+1])/2}) -> mid_pos
floor(mid_pos) -> wid_mid_pos #all wins
rwid_mid_pos_3L=as.matrix(wid_mid_pos[as.numeric(windows)])/10^6       #wins really exsiting in "a"
row.names(rwid_mid_pos_3L)=rownames(a)

#Chr3R
data<-read.table("/home/cmb-12/fs3/hli465/Dpgp/all_samples_Chr3R/pairwise_distance_between_win103_Chr3R_PCAk5",check.names=FALSE)
ind1 <- apply(data, 1, function(x) all(is.na(x)))
a1=data[!ind1,]
ind2 <- apply(a1, 2, function(x) all(is.na(x)))
a2=a1[,!ind2]     
diag(a2)=rep(0,nrow(a2))                                        #delete NAs
a=sqrt(a2)

a <- as.matrix(a)
windows=rownames(a) 
fit2<-cmdscale(a,eig=TRUE, k=2)  
x2 <- fit2$points[,1]
y2 <- fit2$points[,2]
coordinate1<- x2
coordinate1=as.matrix(coordinate1)
row.names(coordinate1)=rownames(a)
coordinate2<- y2
coordinate2=as.matrix(coordinate2)
row.names(coordinate2)=rownames(a)
coordinate1_3R=coordinate1
coordinate2_3R=coordinate2
read.table("/home/cmb-12/fs3/hli465/Dpgp/all_samples_Chr3R/SNP_Pos_Chr3R.txt")[[1]]->snp_pos
floor(length(snp_pos)/1000) -> n
sapply(1:n, function(i){(snp_pos[i*1000]+snp_pos[(i-1)*1000+1])/2}) -> mid_pos
floor(mid_pos) -> wid_mid_pos #all wins
rwid_mid_pos_3R=as.matrix(wid_mid_pos[as.numeric(windows)])/10^6       #wins really exsiting in "a"
row.names(rwid_mid_pos_3R)=rownames(a)

#ChrX
data<-read.table("/home/cmb-12/fs3/hli465/Dpgp/all_samples_ChrX/pairwise_distance_between_win103_ChrX_PCAk5",check.names=FALSE)
ind1 <- apply(data, 1, function(x) all(is.na(x)))
a1=data[!ind1,]
ind2 <- apply(a1, 2, function(x) all(is.na(x)))
a2=a1[,!ind2] 
diag(a2)=rep(0,nrow(a2))                                            #delete NAs
a=sqrt(a2)

a <- as.matrix(a)
windows=rownames(a) 
fit2<-cmdscale(a,eig=TRUE, k=2)  
x2 <- fit2$points[,1]
y2 <- fit2$points[,2]
coordinate1<- x2
coordinate1=as.matrix(coordinate1)
row.names(coordinate1)=rownames(a)
coordinate2<- y2
coordinate2=as.matrix(coordinate2)
row.names(coordinate2)=rownames(a)
coordinate1_X=coordinate1
coordinate2_X=coordinate2
read.table("/home/cmb-12/fs3/hli465/Dpgp/all_samples_ChrX/SNP_Pos_ChrX.txt")[[1]]->snp_pos
floor(length(snp_pos)/1000) -> n
sapply(1:n, function(i){(snp_pos[i*1000]+snp_pos[(i-1)*1000+1])/2}) -> mid_pos
floor(mid_pos) -> wid_mid_pos #all wins
rwid_mid_pos_X=as.matrix(wid_mid_pos[as.numeric(windows)])/10^6       #wins really exsiting in "a"
row.names(rwid_mid_pos_X)=rownames(a)







####put the MDS plots together
pdf(file="/home/cmb-12/fs3/hli465/Dpgp/Fig1_allchr_Together_MDS_plot_PCAk5_win103.pdf",width=6,height=6.5)
layout(matrix(c(1:15),nrow=5,ncol=3,byrow=TRUE),widths=c(1.15,1,1),heights=c(1,1,1,1,1.35))
this.cex <- 0.6
#2L
    
    par(mar=c(1,5.5,0.5,0.5)+.2,xpd=TRUE,mgp=c(2.5,1,0))   
    plot(coordinate1_2L, coordinate2_2L, xlab="",ylab="MDS coordinate 2", pch=20, cex=this.cex)
    mtext("2L",side=2,line=4,font=2)
    par(mar=c(1,3.5,0.5,0.5)+.2,xpd=FALSE,mgp=c(2.5,1,0))
    plot( rwid_mid_pos_2L, coordinate1_2L, xlab="", ylab="MDS coordinate 1", pch=20, cex=this.cex)
    par(mar=c(1,3.5,0.5,0.5)+.2,xpd=FALSE,mgp=c(2.5,1,0))
    abline(v=c(2225744/10^6, 13154180/10^6))
    plot( rwid_mid_pos_2L, coordinate2_2L, xlab="", ylab="MDS coordinate 2",pch=20, cex=this.cex)   
    abline(v=c(2225744/10^6, 13154180/10^6))
 #2R
    
    par(mar=c(1,5.5,0.8,0.5)+.2,xpd=TRUE,mgp=c(2.5,1,0))   
    plot(coordinate1_2R, coordinate2_2R, xlab="",ylab="MDS coordinate 2",pch=20, cex=this.cex)
    mtext("2R",side=2,line=4,font=2)
    par(mar=c(1,3.5,0.8,0.5)+.2,xpd=FALSE,mgp=c(2.5,1,0))
    plot( rwid_mid_pos_2R, coordinate1_2R, xlab="", ylab="MDS coordinate 1", pch=20, cex=this.cex)
    par(mar=c(1,3.5,0.8,0.5)+.2,xpd=FALSE,mgp=c(2.5,1,0))
    abline(v=c(11278659/10^6, 16163839/10^6))
    plot( rwid_mid_pos_2R, coordinate2_2R, xlab="", ylab="MDS coordinate 2",pch=20, cex=this.cex,  )
    abline(v=c(11278659/10^6, 16163839/10^6))
#3L
    
    par(mar=c(1,5.5,0.8,0.5)+.2,xpd=TRUE,mgp=c(2.5,1,0))   
    plot(coordinate1_3L, coordinate2_3L, xlab="",ylab="MDS coordinate 2",pch=20, cex=this.cex)
    mtext("3L",side=2,line=4,font=2)
    par(mar=c(1,3.5,0.8,0.5)+.2,xpd=FALSE,mgp=c(2.5,1,0))
    plot( rwid_mid_pos_3L, coordinate1_3L, xlab="", ylab="MDS coordinate 1", pch=20, cex=this.cex)
    par(mar=c(1,3.5,0.8,0.5)+.2,xpd=FALSE,mgp=c(2.5,1,0))
    abline(v=c(2311085/10^6, 11097829/10^6))
    abline(v=c(3173046/10^6,16301941/10^6),lty=2,col="grey")
    plot( rwid_mid_pos_3L, coordinate2_3L, xlab="", ylab="MDS coordinate 2",pch=20, cex=this.cex)
    abline(v=c(2311085/10^6, 11097829/10^6))
    abline(v=c(3173046/10^6,16301941/10^6),lty=2,col="grey")

#3R
    
    par(mar=c(1,5.5,0.8,0.5)+.2,xpd=TRUE,mgp=c(2.5,1,0))   
    plot(coordinate1_3R, coordinate2_3R, xlab="",ylab="MDS coordinate 2",pch=20, cex=this.cex)
    mtext("3R",side=2,line=4,font=2)
    par(mar=c(1,3.5,0.8,0.5)+.2,xpd=FALSE,mgp=c(2.5,1,0))
    plot( rwid_mid_pos_3R, coordinate1_3R, xlab="", ylab="MDS coordinate 1", pch=20, cex=this.cex)
    par(mar=c(1,3.5,0.8,0.5)+.2,xpd=FALSE,mgp=c(2.5,1,0))
    abline(v=c(7576289/10^6, 21966092/10^6))
    abline(v=c(17232639/10^6,24857019/10^6),lty=2,col="grey")
    abline(v=c(12257931/10^6,20569732/10^6),lty=3,col="grey")
    plot( rwid_mid_pos_3R, coordinate2_3R, xlab="", ylab="MDS coordinate 2",pch=20, cex=this.cex)
    abline(v=c(7576289/10^6, 21966092/10^6))
    abline(v=c(17232639/10^6,24857019/10^6),lty=2,col="grey")
    abline(v=c(12257931/10^6,20569732/10^6),lty=3,col="grey")
#X
    
    par(mar=c(4.1,5.5,0.8,0.5)+.2,xpd=TRUE,mgp=c(2.5,1,0))   
    plot(coordinate1_X, coordinate2_X, xlab="MDS coordinate 1",ylab="MDS coordinate 2",pch=20, cex=this.cex)
    mtext("X",side=2,line=4,font=2)
    par(mar=c(4.1,3.5,0.8,0.5)+.2,xpd=FALSE,mgp=c(2.5,1,0))
    plot( rwid_mid_pos_X, coordinate1_X, xlab="position (Mbp)", ylab="MDS coordinate 1", pch=20, cex=this.cex)
    abline(v=c(13519769/10^6, 19473361/10^6))
    abline(v=c(17722945/10^6,19487744/10^6),lty=2,col="grey")
    par(mar=c(4.1,3.5,0.8,0.5)+.2,xpd=FALSE,mgp=c(2.5,1,0))
    plot( rwid_mid_pos_X, coordinate2_X, xlab="position (Mbp)", ylab="MDS coordinate 2",pch=20, cex=this.cex)
    abline(v=c(13518816/10^6, 19488829/10^6))
    abline(v=c(17722945/10^6,19487744/10^6),lty=2,col="grey")
dev.off()










