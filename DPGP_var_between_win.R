PCAs <- read.table("~/Documents/Drosophila/fluctuation_PCA_win_103_Chr3L.txt",header=TRUE)
len=(nrow(PCAs)-2)/2
k=which(is.na(PCAs[1,]))
PC1s=PCAs[1:len,-k]
a=as.matrix(PC1s)
a=t(a)
A=matrix(0,nrow=nrow(a),ncol=ncol(a))
S=rep(0,nrow(A))
for(i in 1:nrow(A))
{
	if (sum((a[1,]-a[i,])^2)<sum((a[1,]+a[i,])^2))
	{
		S[i]=1
		A[i,]=a[i,]	
	}
	else 
	{
		S[i]=-1
		A[i,]=-a[i,]	
	}
}
varfunction<-function(data){
	var=sum((data-mean(data))^2)/(nrow(A))
	return(var)
}
b=apply(A,2,varfunction)
mb=mean(b)
