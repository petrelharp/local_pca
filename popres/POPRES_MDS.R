a<-read.table("quick_method_pairwise_distance_between_win_100_chr1")
k=which(is.na(a),arr.ind=TRUE)
a[k]=0
a=a+t(a)
a <- as.matrix(a)
a=sqrt(a)
fit2d<-cmdscale(a,eig=TRUE, k=2)  # MDS_2D#
x <- fit2d$points[,1]
y <- fit2d$points[,2]
pdf(file="MDS_2D_win100_chr1.pdf")
layout(1:2)
par(mar=c(4,3,1,1)+0.1)
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", main="MDS_2D_chr1", col=rainbow(2*nrow(a)))
plot(1:nrow(a),col=rainbow(2*nrow(a)))
dev.off()

fit1d<-cmdscale(a,eig=TRUE, k=1)  # MDS_1D#
x <- fit1d$points[,1]
pdf(file="MDS_1D_win100_chr1.pdf")
plot(x,ylab="Coordinate 1", main="MDS_1D_chr2")
dev.off()
