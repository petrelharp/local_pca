
util.plot_medicago <- function(chr) {

original.wd <- getwd()
setwd("/home/cmb-12/fs3/hli465/FormatVCFbgz")

data<-read.table(gettextf("pairwise_distance_between_win104_%s_PCAk2",chr),check.names=FALSE)
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

#MDS_SNP_Pos (we choosed this one~)
read.table(gettextf("SNP_pos_%s", chr))[[1]]->snp_pos
floor(length(snp_pos)/10^4) -> n
sapply(1:n, function(i){(snp_pos[i*10^4]+snp_pos[(i-1)*10^4+1])/2}) -> mid_pos
floor(mid_pos) -> wid_mid_pos #all wins
rwid_mid_pos=as.matrix(wid_mid_pos[as.numeric(windows)])/10^6     #wins really exsiting in "a"
row.names(rwid_mid_pos)=rownames(a)
missings<-read.table(gettextf("missings_for_each_pos_%s.txt",chr))
mis_win=sapply(1:n, function(i){mean(missings[((i-1)*10000+1):(i*10000),])})
mis_win=mis_win[as.numeric(windows)]

result <- list(coordinate1=coordinate1,coordinate2=coordinate2,rwid_mid_pos=rwid_mid_pos,mis_win=mis_win)
return(result)

setwd(original.wd)

}


rchr1=util.plot_medicago("chr1")     #result of chr1
rchr2=util.plot_medicago("chr2")
rchr3=util.plot_medicago("chr3")
rchr4=util.plot_medicago("chr4")     
rchr5=util.plot_medicago("chr5")     
rchr6=util.plot_medicago("chr6")    
rchr7=util.plot_medicago("chr7")    
rchr8=util.plot_medicago("chr8")     

pdf(file="/home/cmb-12/fs3/hli465/FormatVCFbgz/Medicago_missing_against_MDS1_win104.pdf",width=5,height=7)
layout(matrix(c(1:16), nrow=8,ncol=2,byrow=TRUE),widths=c(1.1,1),heights=c(rep(1,7),1.3))
this.cex <- 0.5
#chr1
par(mar=c(1.5,4.5,0.6,0.5),xpd=TRUE)
plot(rchr1$rwid_mid_pos, rchr1$mis_win, xlab="", ylab="missingness",pch=20,cex=this.cex,cex.lab=0.65,cex.axis=0.65,mgp=c(2,1,0),mgp=c(2,1,0))
mtext("1",side=2,line=3,font=1.5)
par(mar=c(1.5,3,0.6,0.5),xpd=TRUE)
plot(rchr1$mis_win,rchr1$coordinate1,xlab="",ylab="MDS coordinate1",cex=this.cex,cex.lab=0.65,cex.axis=0.65,mgp=c(2,1,0))

#chr2
par(mar=c(1.5,4.5,0.6,0.5),xpd=TRUE)
plot(rchr2$rwid_mid_pos, rchr2$mis_win, xlab="", ylab="missingness", pch=20,cex=this.cex,cex.lab=0.65,cex.axis=0.65,mgp=c(2,1,0))
mtext("2",side=2,line=3,font=1.5)
par(mar=c(1.5,3,0.6,0.5),xpd=TRUE)
plot(rchr2$mis_win,rchr2$coordinate1,xlab="",ylab="MDS coordinate1",,cex=this.cex,cex.lab=0.65,cex.axis=0.65,mgp=c(2,1,0))

#chr3
par(mar=c(1.5,4.5,0.6,0.5),xpd=TRUE)
plot(rchr3$rwid_mid_pos, rchr3$mis_win, xlab="", ylab="missingness", pch=20,cex=this.cex,cex.lab=0.65,cex.axis=0.65,mgp=c(2,1,0))
mtext("3",side=2,line=3,font=1.5)
par(mar=c(1.5,3,0.6,0.5),xpd=TRUE)
plot(rchr3$mis_win,rchr3$coordinate1,xlab="",ylab="MDS coordinate1",,cex=this.cex,cex.lab=0.65,cex.axis=0.65,mgp=c(2,1,0))

#chr4
par(mar=c(1.5,4.5,0.6,0.5),xpd=TRUE)
plot(rchr4$rwid_mid_pos, rchr4$mis_win, xlab="", ylab="missingness", pch=20,cex=this.cex,cex.lab=0.65,cex.axis=0.65,mgp=c(2,1,0))
mtext("4",side=2,line=3,font=1.5)
par(mar=c(1.5,3,0.6,0.5),xpd=TRUE)
plot(rchr4$mis_win,rchr4$coordinate1,xlab="",ylab="MDS coordinate1",,cex=this.cex,cex.lab=0.65,cex.axis=0.65,mgp=c(2,1,0))

#chr5
par(mar=c(1.5,4.5,0.6,0.5),xpd=TRUE)
plot(rchr5$rwid_mid_pos, rchr5$mis_win, xlab="", ylab="missingness", pch=20,cex=this.cex,cex.lab=0.65,cex.axis=0.65,mgp=c(2,1,0))
mtext("5",side=2,line=3,font=1.5)
par(mar=c(1.5,3,0.6,0.5),xpd=TRUE)
plot(rchr5$mis_win,rchr5$coordinate1,xlab="",ylab="MDS coordinate1",,cex=this.cex,cex.lab=0.65,cex.axis=0.65,mgp=c(2,1,0))

#chr6
par(mar=c(1.5,4.5,0.6,0.5),xpd=TRUE)
plot(rchr6$rwid_mid_pos, rchr6$mis_win, xlab="", ylab="missingness", pch=20,cex=this.cex,cex.lab=0.65,cex.axis=0.65,mgp=c(2,1,0))
mtext("6",side=2,line=3,font=1.5)
par(mar=c(1.5,3,0.6,0.5),xpd=TRUE)
plot(rchr6$mis_win,rchr6$coordinate1,xlab="",ylab="MDS coordinate1",,cex=this.cex,cex.lab=0.65,cex.axis=0.65,mgp=c(2,1,0))

#chr7
par(mar=c(1.5,4.5,0.6,0.5),xpd=TRUE)
plot(rchr7$rwid_mid_pos, rchr7$mis_win, xlab="", ylab="missingness", pch=20,cex=this.cex,cex.lab=0.65,cex.axis=0.65,mgp=c(2,1,0))
mtext("7",side=2,line=3,font=1.5)
par(mar=c(1.5,3,0.6,0.5),xpd=TRUE)
plot(rchr7$mis_win,rchr7$coordinate1,xlab="",ylab="MDS coordinate1",cex=this.cex,cex.lab=0.65,cex.axis=0.65,mgp=c(2,1,0))

#chr8
par(mar=c(3,4.5,0.6,0.5),xpd=TRUE)
plot(rchr8$rwid_mid_pos, rchr8$mis_win, xlab="position (Mbp)", ylab="missingness", pch=20,cex=this.cex,cex.lab=0.65,cex.axis=0.65,mgp=c(2,1,0))
mtext("8",side=2,line=3,font=1.5)
par(mar=c(3,3,0.6,0.5),xpd=TRUE)
plot(rchr8$mis_win,rchr8$coordinate1,xlab="missingness",ylab="MDS coordinate1",,cex=this.cex,cex.lab=0.65,cex.axis=0.65,mgp=c(2,1,0))


dev.off()





