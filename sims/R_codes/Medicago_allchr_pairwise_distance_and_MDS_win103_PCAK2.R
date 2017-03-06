#win103 PCAK2
fn=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8")
win2 <- 10^3
for(i in 1:length(fn))
{
    data <- read.table(paste("/home/cmb-12/fs3/hli465/FormatVCFbgz/",fn[i],sep=""),na.strings=".",stringsAsFactors=FALSE)
    data<-data[,-(1:2)]
    n<-ncol(data)/2
    coded<-data[,2*(1:n)]+data[,2*(1:n)-1]
    eigenstuff <- eigen_windows(coded, win=10^3, k=2,recode=FALSE)
    windist <- pc_dist(eigenstuff, npc=2)
    write.table(windist,file=paste("/home/cmb-12/fs3/hli465/FormatVCFbgz/pairwise_distance_between_win103",fn[i],"PCAk2",sep="_"),sep="\t")
    rm(data,coded,eigenstuff,windist)
    gc()
}

util.plot_medicago <- function(chr) {


setwd("/home/cmb-12/fs3/hli465/FormatVCFbgz")

data<-read.table(gettextf("pairwise_distance_between_win103_%s_PCAk2",chr),check.names=FALSE)
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

read.table(gettextf("SNP_pos_%s", chr))[[1]]->snp_pos
floor(length(snp_pos)/10^3) -> n
sapply(1:n, function(i){(snp_pos[i*10^3]+snp_pos[(i-1)*10^3+1])/2}) -> mid_pos
floor(mid_pos) -> wid_mid_pos #all wins
rwid_mid_pos=as.matrix(wid_mid_pos[as.numeric(windows)])/10^6     #wins really exsiting in "a"
row.names(rwid_mid_pos)=rownames(a)

result <- list(coordinate1=coordinate1,coordinate2=coordinate2,rwid_mid_pos=rwid_mid_pos)
return(result)

}


rchr1=util.plot_medicago("chr1")     #result of chr1
rchr2=util.plot_medicago("chr2")
rchr3=util.plot_medicago("chr3")
rchr4=util.plot_medicago("chr4")     
rchr5=util.plot_medicago("chr5")     
rchr6=util.plot_medicago("chr6")    
rchr7=util.plot_medicago("chr7")    
rchr8=util.plot_medicago("chr8")     

pdf(file="/home/cmb-12/fs3/hli465/FormatVCFbgz/Medicago_MDS_plot_allchr_win103_PCAk2.pdf",width=5,height=7)
layout(matrix(c(1:24), nrow=8,ncol=3,byrow=TRUE),widths=c(1.1,1,1),heights=c(rep(1,7),1.3))
this.cex <- 0.5
#chr1
par(mar=c(1.5,4.5,0.6,0.5),xpd=TRUE)
plot(rchr1$coordinate1, rchr1$coordinate2, xlab="", ylab="MDS coordinate2",pch=20,cex=this.cex,cex.lab=0.65,cex.axis=0.65,mgp=c(2,1,0),mgp=c(2,1,0))
mtext("1",side=2,line=3,font=1.5)
par(mar=c(1.5,3,0.6,0.5),xpd=TRUE)
plot(rchr1$rwid_mid_pos,rchr1$coordinate1,xlab="",ylab="MDS coordinate1",cex=this.cex,cex.lab=0.65,cex.axis=0.65,mgp=c(2,1,0))
par(mar=c(1.5,3,0.6,0.5),xpd=TRUE)
plot(rchr1$rwid_mid_pos, rchr1$coordinate2,xlab="", ylab="MDS coordinate2",cex=this.cex,cex.lab=0.65,cex.axis=0.65,mgp=c(2,1,0))

#chr2
par(mar=c(1.5,4.5,0.6,0.5),xpd=TRUE)
plot(rchr2$coordinate1, rchr2$coordinate2, xlab="", ylab="MDS coordinate2", pch=20,cex=this.cex,cex.lab=0.65,cex.axis=0.65,mgp=c(2,1,0))
mtext("2",side=2,line=3,font=1.5)
par(mar=c(1.5,3,0.6,0.5),xpd=TRUE)
plot(rchr2$rwid_mid_pos,rchr2$coordinate1,xlab="",ylab="MDS coordinate1",,cex=this.cex,cex.lab=0.65,cex.axis=0.65,mgp=c(2,1,0))
par(mar=c(1.5,3,0.6,0.5),xpd=TRUE)
plot(rchr2$rwid_mid_pos,rchr2$coordinate2,xlab="", ylab="MDS coordinate2",cex=this.cex,cex.lab=0.65,cex.axis=0.65,mgp=c(2,1,0))

#chr3
par(mar=c(1.5,4.5,0.6,0.5),xpd=TRUE)
plot(rchr3$coordinate1, rchr3$coordinate2, xlab="", ylab="MDS coordinate2", pch=20,cex=this.cex,cex.lab=0.65,cex.axis=0.65,mgp=c(2,1,0))
mtext("3",side=2,line=3,font=1.5)
par(mar=c(1.5,3,0.6,0.5),xpd=TRUE)
plot(rchr3$rwid_mid_pos,rchr3$coordinate1,xlab="",ylab="MDS coordinate1",,cex=this.cex,cex.lab=0.65,cex.axis=0.65,mgp=c(2,1,0))
par(mar=c(1.5,3,0.6,0.5),xpd=TRUE)
plot(rchr3$rwid_mid_pos,rchr3$coordinate2,xlab="", ylab="MDS coordinate2",cex=this.cex,cex.lab=0.65,cex.axis=0.65,mgp=c(2,1,0))

#chr4
par(mar=c(1.5,4.5,0.6,0.5),xpd=TRUE)
plot(rchr4$coordinate1, rchr4$coordinate2, xlab="", ylab="MDS coordinate2", pch=20,cex=this.cex,cex.lab=0.65,cex.axis=0.65,mgp=c(2,1,0))
mtext("4",side=2,line=3,font=1.5)
par(mar=c(1.5,3,0.6,0.5),xpd=TRUE)
plot(rchr4$rwid_mid_pos,rchr4$coordinate1,xlab="",ylab="MDS coordinate1",,cex=this.cex,cex.lab=0.65,cex.axis=0.65,mgp=c(2,1,0))
par(mar=c(1.5,3,0.6,0.5),xpd=TRUE)
plot(rchr4$rwid_mid_pos,rchr4$coordinate2,xlab="", ylab="MDS coordinate2",cex=this.cex,cex.lab=0.65,cex.axis=0.65,mgp=c(2,1,0))

#chr5
par(mar=c(1.5,4.5,0.6,0.5),xpd=TRUE)
plot(rchr5$coordinate1, rchr5$coordinate2, xlab="", ylab="MDS coordinate2", pch=20,cex=this.cex,cex.lab=0.65,cex.axis=0.65,mgp=c(2,1,0))
mtext("5",side=2,line=3,font=1.5)
par(mar=c(1.5,3,0.6,0.5),xpd=TRUE)
plot(rchr5$rwid_mid_pos,rchr5$coordinate1,xlab="",ylab="MDS coordinate1",,cex=this.cex,cex.lab=0.65,cex.axis=0.65,mgp=c(2,1,0))
par(mar=c(1.5,3,0.6,0.5),xpd=TRUE)
plot(rchr5$rwid_mid_pos,rchr5$coordinate2,xlab="", ylab="MDS coordinate2",cex=this.cex,cex.lab=0.65,cex.axis=0.65,mgp=c(2,1,0))

#chr6
par(mar=c(1.5,4.5,0.6,0.5),xpd=TRUE)
plot(rchr6$coordinate1, rchr6$coordinate2, xlab="", ylab="MDS coordinate2", pch=20,cex=this.cex,cex.lab=0.65,cex.axis=0.65,mgp=c(2,1,0))
mtext("6",side=2,line=3,font=1.5)
par(mar=c(1.5,3,0.6,0.5),xpd=TRUE)
plot(rchr6$rwid_mid_pos,rchr6$coordinate1,xlab="",ylab="MDS coordinate1",,cex=this.cex,cex.lab=0.65,cex.axis=0.65,mgp=c(2,1,0))
par(mar=c(1.5,3,0.6,0.5),xpd=TRUE)
plot(rchr6$rwid_mid_pos,rchr6$coordinate2,xlab="", ylab="MDS coordinate2",cex=this.cex,cex.lab=0.65,cex.axis=0.65,mgp=c(2,1,0))

#chr7
par(mar=c(1.5,4.5,0.6,0.5),xpd=TRUE)
plot(rchr7$coordinate1, rchr7$coordinate2, xlab="", ylab="MDS coordinate2", pch=20,cex=this.cex,cex.lab=0.65,cex.axis=0.65,mgp=c(2,1,0))
mtext("7",side=2,line=3,font=1.5)
par(mar=c(1.5,3,0.6,0.5),xpd=TRUE)
plot(rchr7$rwid_mid_pos,rchr7$coordinate1,xlab="",ylab="MDS coordinate1",cex=this.cex,cex.lab=0.65,cex.axis=0.65,mgp=c(2,1,0))
par(mar=c(1.5,3,0.6,0.5),xpd=TRUE)
plot(rchr7$rwid_mid_pos,rchr7$coordinate2,xlab="", ylab="MDS coordinate2",cex=this.cex,cex.lab=0.65,cex.axis=0.65,mgp=c(2,1,0))

#chr8
par(mar=c(3,4.5,0.6,0.5),xpd=TRUE)
plot(rchr8$coordinate1, rchr8$coordinate2, xlab="MDS coordinate1", ylab="MDS coordinate2", pch=20,cex=this.cex,cex.lab=0.65,cex.axis=0.65,mgp=c(2,1,0))
mtext("8",side=2,line=3,font=1.5)
par(mar=c(3,3,0.6,0.5),xpd=TRUE)
plot(rchr8$rwid_mid_pos,rchr8$coordinate1,xlab="position (Mbp)",ylab="MDS coordinate1",,cex=this.cex,cex.lab=0.65,cex.axis=0.65,mgp=c(2,1,0))
par(mar=c(3,3,0.6,0.5),xpd=TRUE)
plot(rchr8$rwid_mid_pos,rchr8$coordinate2,xlab="position (Mbp)", ylab="MDS coordinate2",cex=this.cex,cex.lab=0.65,cex.axis=0.65,mgp=c(2,1,0))


dev.off()





