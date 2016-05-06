a<-read.table("~/Documents/Drosophila/Chr3L/pairwise_distance_between_win_103_Chr3L")
a <- as.matrix(a)
count=floor(nrow(a)*0.05)
fit1<-cmdscale(a,eig=TRUE, k=2)  
x <- fit1$points[,1]
y <- fit1$points[,2]
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", main="MDS_2D", col=rainbow(2*nrow(a)))
b=identify(x,y,n=3,labels=row.names(a))
b
get.neibor<-function(x,y){
	     temp=sort(y[,x])
	     neibor=as.numeric(names(temp)[1:count])
	     return(neibor)
}
neibors<-sapply(b[1:3],get.neibor,y=a)
write.table(neibors,file="~/Documents/Drosophila/Chr3L/neibors_win103_Chr3L.txt",sep="\t")
