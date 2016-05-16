setwd("~/Documents/Drosophila/Chr2L")
a<-read.table("pairwise_distance_between_win_103_Chr2L",check.names=FALSE)
a <- as.matrix(a)
b1<-as.matrix(read.table("neibors_win103_Chr2L_ordered.txt"))
windows=rownames(a) 
fit2<-cmdscale(a,eig=TRUE, k=2)  
x2 <- fit2$points[,1]
y2 <- fit2$points[,2]
temp1=as.matrix(cbind(x2,y2))
dist1=as.matrix(dist(temp1))
rownames(dist1)=rownames(a)
colnames(dist1)=colnames(a)
count=floor(nrow(a)*0.05)
get.neibor<-function(x,y){
	     temp=sort(y[,x])
	     neibor=as.numeric(names(temp)[1:count])
	     return(neibor)
}
neibors<-sapply(as.character(b1[1,]),get.neibor,y=dist1)
write.table(neibors,file="neibors_win103_Chr2L_byMDS_ordered.txt",sep="\t")

library(RColorBrewer)
display.brewer.all(3,colorblindFriendly=TRUE)
cols <- brewer.pal(3,"Dark2")
b<-neibors
coordinate1<- x2
coordinate1=as.matrix(coordinate1)
row.names(coordinate1)=rownames(a)
coordinate2<- y2
coordinate2=as.matrix(coordinate2)
row.names(coordinate2)=rownames(a)
color=rep("black",nrow(coordinate1))
Pch=rep(1,nrow(coordinate1))
color[rownames(coordinate1)%in%as.character(b[,1])]=cols[1]
Pch[rownames(coordinate1)%in%as.character(b[,1])]=15
color[rownames(coordinate1)%in%as.character(b[,2])]=cols[2]
Pch[rownames(coordinate1)%in%as.character(b[,2])]=17
color[rownames(coordinate1)%in%as.character(b[,3])]=cols[3]
Pch[rownames(coordinate1)%in%as.character(b[,3])]=19

#MDS_winID
#pdf(file="Fig1_Together_MDS_plot_Chr2L_try3.pdf",width=8,height=8)
#layout(matrix(c(1,1,2,3), nrow=2,ncol=2,byrow=TRUE),widths=c(2,2,2,2),heights=c(1,1,2,2))
# par(mar=c(4,4,4,2),xpd=TRUE)
# plot(x2, y2, xlab="Coordinate1", ylab="Coordinate2", main="MDS_2D_Chr2L", col=color,pch=Pch)

# plot(windows,coordinate1,xlab="windows_Chr2L", main="MDS_1_Chr2L")
# points(c(b[,2]),coordinate1[as.character(c(b[,2])),1],col=cols[3],pch=15)
# points(c(b[,1]),coordinate1[as.character(c(b[,1])),1],col=cols[1],pch=17)
# points(c(b[,3]),coordinate1[as.character(c(b[,3])),1],col=cols[2],pch=19)

# plot(windows,coordinate2,xlab="windows_Chr2L", main="MDS_2_Chr2L")
# points(c(b[,2]),coordinate2[as.character(c(b[,2])),1],col=cols[3],pch=15)
# points(c(b[,1]),coordinate2[as.character(c(b[,1])),1],col=cols[1],pch=17)
# points(c(b[,3]),coordinate2[as.character(c(b[,3])),1],col=cols[2],pch=19)
# dev.off()



#MDS_SNP_Pos (we choosed this one~)
read.table("SNP_Pos_Chr2L.txt")[[1]]->snp_pos
floor(length(snp_pos)/1000) -> n
sapply(1:n, function(i){(snp_pos[i*1000]+snp_pos[(i-1)*1000+1])/2}) -> mid_pos
floor(mid_pos) -> wid_mid_pos #all wins
rwid_mid_pos=as.matrix(wid_mid_pos[as.numeric(windows)])       #wins really exsiting in "a"
row.names(rwid_mid_pos)=rownames(a)

pdf(file="Fig1_Together_MDS_plot_Chr2L_final_abline.pdf",width=10,height=3)
layout(matrix(c(1,2,3), nrow=1,ncol=3,byrow=TRUE))
par(mar=c(4.5,4,4,2),xpd=TRUE)
plot(x2, y2, xlab="coordinate1", ylab="coordinate2", col=color,pch=Pch)

plot(rwid_mid_pos,coordinate1,xlab="position on Chr2L (bp)",main="Drosophila Chromosome 2L")
points(rwid_mid_pos[as.character(b[,1]),1],coordinate1[as.character(b[,1]),1],col=cols[1],pch=15)
points(rwid_mid_pos[as.character(b[,2]),1],coordinate1[as.character(b[,2]),1],col=cols[2],pch=17)
points(rwid_mid_pos[as.character(b[,3]),1],coordinate1[as.character(b[,3]),1],col=cols[3],pch=19)
abline(v=c(2225744, 13154180))

plot(rwid_mid_pos,coordinate2,xlab="position on Chr2L (bp)")
points(rwid_mid_pos[as.character(b[,1]),1],coordinate2[as.character(b[,1]),1],col=cols[1],pch=15)
points(rwid_mid_pos[as.character(b[,2]),1],coordinate2[as.character(b[,2]),1],col=cols[2],pch=17)
points(rwid_mid_pos[as.character(b[,3]),1],coordinate2[as.character(b[,3]),1],col=cols[3],pch=19)
abline(v=c(2225744, 13154180))
dev.off()





#MDS_SNP_pos no cover points
pdf(file="Fig1_Together_MDS_plot_Chr2L_final1.pdf",width=10,height=3)
layout(matrix(c(1,2,3), nrow=1,ncol=3,byrow=TRUE))
par(mar=c(4.5,4,4,2),xpd=TRUE)
plot(x2, y2, xlab="coordinate1", ylab="coordinate2", col=color,pch=Pch)
plot(rwid_mid_pos,coordinate1,xlab="position on chromosome (bp)",main="Drosophila Chromosome 2R",col=color,pch=Pch)
plot(rwid_mid_pos,coordinate2,xlab="position on chromosome (bp)",col=color,pch=Pch)
dev.off()











