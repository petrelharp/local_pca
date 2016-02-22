setwd("~/Documents/Drosophila/Chr2L")
a<-read.table("pairwise_distance_between_win_103_Chr2L")
a <- as.matrix(a)
b<-read.table("neibors_win103_Chr2L.txt")
fit1<-cmdscale(a,eig=TRUE, k=1)
windows=rownames(a) 
windows=as.numeric(windows) 
coordinate1<- fit1$points[,1]
fit2<-cmdscale(a,eig=TRUE, k=2)  
x2 <- fit2$points[,1]
y2 <- fit2$points[,2]
coordinate1=as.matrix(coordinate1)
row.names(coordinate1)=rownames(a)
pdf(file="Fig1_Together_MDS_plot_Chr2L.pdf",width=12,height=6)
layout(matrix(c(1,2,3,3), nrow=2),heights=c(2,1,3))
plot(x2, y2, xlab="Coordinate1", ylab="Coordinate2", main="MDS_2D_Chr2L", col=rainbow(2*nrow(a)))
plot(1:nrow(a),y=rep(1,nrow(a)),col=rainbow(2*nrow(a)),xlab=' ',ylab=' ',yaxt='n')
plot(windows,coordinate1,xlab="windows_chr2L", main="MDS_1D_Chr2L")
points(c(b[,3]),coordinate1[as.character(c(b[,3])),1],col="blue",pch=15)
points(c(b[,2]),coordinate1[as.character(c(b[,2])),1],col="green",pch=17)
points(c(b[,1]),coordinate1[as.character(c(b[,1])),1],col="red",pch=19)
dev.off()


