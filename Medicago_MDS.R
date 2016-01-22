setwd("~/Documents/Medicago_Project/gene_density_related")
a<-read.table("quick_method_pairwise_distance_between_win_104_all_chr_medicago")
k=which(is.na(a),arr.ind=TRUE)
a[k]=0
a=a+t(a)
a <- as.matrix(a)
a=sqrt(a)

fit1<-cmdscale(a,eig=TRUE, k=1)  # MDS_1D#
x <- fit1$points[,1]
pdf(file="MDS_1D_win104_all_chr.pdf")
layout(1:2)
par(mar=c(4,3,1,1)+0.1)
plot(x,ylab="Coordinate 1", main="MDS_1D_all_chr", col=rainbow(2*nrow(a)))
plot(1:nrow(a),col=rainbow(2*nrow(a)))
dev.off()


#col MDS1 by chr
read.table("SNP_pos_chr1")[[1]]->snp_pos
floor(length(snp_pos)/10000) -> n1
read.table("SNP_pos_chr2")[[1]]->snp_pos
floor(length(snp_pos)/10000) -> n2
read.table("SNP_pos_chr3")[[1]]->snp_pos
floor(length(snp_pos)/10000) -> n3
read.table("SNP_pos_chr4")[[1]]->snp_pos
floor(length(snp_pos)/10000) -> n4
read.table("SNP_pos_chr5")[[1]]->snp_pos
floor(length(snp_pos)/10000) -> n5
read.table("SNP_pos_chr6")[[1]]->snp_pos
floor(length(snp_pos)/10000) -> n6
read.table("SNP_pos_chr7")[[1]]->snp_pos
floor(length(snp_pos)/10000) -> n7
read.table("SNP_pos_chr8")[[1]]->snp_pos
floor(length(snp_pos)/10000) -> n8
w=c(n1,n2,n3,n4,n5,n6,n7,n8)
w1=w
w2=rep(1,w[1])
for(i in 2:8) {wi=rep(i,w[i]);w2=c(w2,wi)}

pdf(file="MDS_1D_win104_col_by_chr_medicago.pdf",width=10,height=4)
plot(x,ylab="Coordinate 1", main="MDS_1D_all_chr", pch=w2, col=rainbow(8)[w2],cex=0.5)
legend(4000,0.4,pch=1:8,col=rainbow(8),legend=levels(factor(w2)))
dev.off()

# MDS_for_each_chr_from_all_chr_pairwise_distance

pdf(file="MDS_1D_win104_chr1_based_all_chr.pdf")
plot(x[1:w[1]],ylab="Coordinate 1", main="MDS_1D_chr1", pch=1, col=rainbow(8)[1])
dev.off()

pdf(file="MDS_1D_win104_chr2_based_all_chr.pdf")
plot(x[(1+w[1]):sum(w[1:2])],ylab="Coordinate 1", main="MDS_1D_chr2", pch=2, col=rainbow(8)[2])
dev.off()

pdf(file="MDS_1D_win104_chr3_based_all_chr.pdf")
plot(x[sum(w[1:2],1):sum(w[1:3])],ylab="Coordinate 1", main="MDS_1D_chr3", pch=3, col=rainbow(8)[3])
dev.off()

pdf(file="MDS_1D_win104_chr4_based_all_chr.pdf")
plot(x[sum(w[1:3],1):sum(w[1:4])],ylab="Coordinate 1", main="MDS_1D_chr4", pch=4, col=rainbow(8)[4])
dev.off()

pdf(file="MDS_1D_win104_chr5_based_all_chr.pdf")
plot(x[sum(w[1:4],1):sum(w[1:5])],ylab="Coordinate 1", main="MDS_1D_chr5", pch=5, col=rainbow(8)[5])
dev.off()

pdf(file="MDS_1D_win104_chr6_based_all_chr.pdf")
plot(x[sum(w[1:5],1):sum(w[1:6])],ylab="Coordinate 1", main="MDS_1D_chr6", pch=6, col=rainbow(8)[6])
dev.off()

pdf(file="MDS_1D_win104_chr7_based_all_chr.pdf")
plot(x[sum(w[1:6],1):sum(w[1:7])],ylab="Coordinate 1", main="MDS_1D_chr7", pch=7, col=rainbow(8)[7])
dev.off()

pdf(file="MDS_1D_win104_chr8_based_all_chr.pdf")
plot(x[sum(w[1:7],1):sum(w[1:8])],ylab="Coordinate 1", main="MDS_1D_chr8", pch=8, col=rainbow(8)[8])
dev.off()

#MDS_for_each_chr_from_all_chr_pairwise_distance_in_bp

read.table("SNP_pos_chr1")[[1]]->snp_pos
floor(length(snp_pos)/10000) -> n
sapply(1:n, function(i){(snp_pos[i*10000]+snp_pos[(i-1)*10000+1])/2}) -> mid_pos
floor(mid_pos) -> wid_mid_pos
pdf(file="MDS_1D_win104_chr1_based_all_chr_SNP_pos.pdf")
plot(wid_mid_pos,x[1:w[1]],ylab="Coordinate 1", main="MDS_1D_chr1", pch=1, col=rainbow(8)[1])
dev.off()


read.table("SNP_pos_chr2")[[1]]->snp_pos
floor(length(snp_pos)/10000) -> n
sapply(1:n, function(i){(snp_pos[i*10000]+snp_pos[(i-1)*10000+1])/2}) -> mid_pos
floor(mid_pos) -> wid_mid_pos
pdf(file="MDS_1D_win104_chr2_based_all_chr_SNP_pos.pdf")
plot(wid_mid_pos,x[(1+w[1]):sum(w[1:2])],ylab="Coordinate 1", main="MDS_1D_chr2", pch=2, col=rainbow(8)[2])
dev.off()


read.table("SNP_pos_chr3")[[1]]->snp_pos
floor(length(snp_pos)/10000) -> n
sapply(1:n, function(i){(snp_pos[i*10000]+snp_pos[(i-1)*10000+1])/2}) -> mid_pos
floor(mid_pos) -> wid_mid_pos
pdf(file="MDS_1D_win104_chr3_based_all_chr_SNP_pos.pdf")
plot(wid_mid_pos,x[sum(w[1:2],1):sum(w[1:3])],ylab="Coordinate 1", main="MDS_1D_chr3", pch=3, col=rainbow(8)[3])
dev.off()

read.table("SNP_pos_chr4")[[1]]->snp_pos
floor(length(snp_pos)/10000) -> n
sapply(1:n, function(i){(snp_pos[i*10000]+snp_pos[(i-1)*10000+1])/2}) -> mid_pos
floor(mid_pos) -> wid_mid_pos
pdf(file="MDS_1D_win104_chr4_based_all_chr_SNP_pos.pdf")
plot(wid_mid_pos,x[sum(w[1:3],1):sum(w[1:4])],ylab="Coordinate 1", main="MDS_1D_chr4", pch=4, col=rainbow(8)[4])
dev.off()

read.table("SNP_pos_chr5")[[1]]->snp_pos
floor(length(snp_pos)/10000) -> n
sapply(1:n, function(i){(snp_pos[i*10000]+snp_pos[(i-1)*10000+1])/2}) -> mid_pos
floor(mid_pos) -> wid_mid_pos
pdf(file="MDS_1D_win104_chr5_based_all_chr_SNP_pos.pdf")
plot(wid_mid_pos,x[sum(w[1:4],1):sum(w[1:5])],ylab="Coordinate 1", main="MDS_1D_chr5", pch=5, col=rainbow(8)[5])
dev.off()

read.table("SNP_pos_chr6")[[1]]->snp_pos
floor(length(snp_pos)/10000) -> n
sapply(1:n, function(i){(snp_pos[i*10000]+snp_pos[(i-1)*10000+1])/2}) -> mid_pos
floor(mid_pos) -> wid_mid_pos
pdf(file="MDS_1D_win104_chr6_based_all_chr_SNP_pos.pdf")
plot(mid_pos,x[sum(w[1:5],1):sum(w[1:6])],ylab="Coordinate 1", main="MDS_1D_chr6", pch=6, col=rainbow(8)[6])
dev.off()

read.table("SNP_pos_chr7")[[1]]->snp_pos
floor(length(snp_pos)/10000) -> n
sapply(1:n, function(i){(snp_pos[i*10000]+snp_pos[(i-1)*10000+1])/2}) -> mid_pos
floor(mid_pos) -> wid_mid_pos
pdf(file="MDS_1D_win104_chr7_based_all_chr_SNP_pos.pdf")
plot(wid_mid_pos,x[sum(w[1:6],1):sum(w[1:7])],ylab="Coordinate 1", main="MDS_1D_chr7", pch=7, col=rainbow(8)[7])
dev.off()

read.table("SNP_pos_chr8")[[1]]->snp_pos
floor(length(snp_pos)/10000) -> n
sapply(1:n, function(i){(snp_pos[i*10000]+snp_pos[(i-1)*10000+1])/2}) -> mid_pos
floor(mid_pos) -> wid_mid_pos
pdf(file="MDS_1D_win104_chr8_based_all_chr_SNP_pos.pdf")
plot(wid_mid_pos,x[sum(w[1:7],1):sum(w[1:8])],ylab="Coordinate 1", main="MDS_1D_chr8", pch=8, col=rainbow(8)[8])
dev.off()



#MDS_against_gene_count_win104

setwd("~/Documents/Medicago_Project/gene_density_related")
a<-read.table("quick_method_pairwise_distance_between_win_104_all_chr_medicago")
k=which(is.na(a),arr.ind=TRUE)
a[k]=0
a=a+t(a)
a <- as.matrix(a)
a=sqrt(a)
fit1<-cmdscale(a,eig=TRUE, k=1)  # MDS_1D#
x <- fit1$points[,1]
read.table("SNP_pos_chr1")[[1]]->snp_pos
floor(length(snp_pos)/10000) -> n1
read.table("SNP_pos_chr2")[[1]]->snp_pos
floor(length(snp_pos)/10000) -> n2
read.table("SNP_pos_chr3")[[1]]->snp_pos
floor(length(snp_pos)/10000) -> n3
read.table("SNP_pos_chr4")[[1]]->snp_pos
floor(length(snp_pos)/10000) -> n4
read.table("SNP_pos_chr5")[[1]]->snp_pos
floor(length(snp_pos)/10000) -> n5
read.table("SNP_pos_chr6")[[1]]->snp_pos
floor(length(snp_pos)/10000) -> n6
read.table("SNP_pos_chr7")[[1]]->snp_pos
floor(length(snp_pos)/10000) -> n7
read.table("SNP_pos_chr8")[[1]]->snp_pos
floor(length(snp_pos)/10000) -> n8
w=c(n1,n2,n3,n4,n5,n6,n7,n8)
gcount=read.table("~/Documents/Medicago_Project/gene_density_related/gene_and_TE_based_gene/chr1_TE_based_gene_starts_count_win104")
gene_count=gcount[-1,3]
gene_count=as.vector(gene_count[-length(gene_count)])
pdf(file="MDS_1D_against_TE_based_gene_count_chr1_based_all_chr_win104.pdf")
plot(gene_count,x[1:w[1]],ylab="Coordinate 1", main="MDS_1D_against_TE_based_gene_count_chr1_win104", pch=1, col=rainbow(8)[1])
dev.off()

gcount=read.table("~/Documents/Medicago_Project/gene_density_related/gene_and_TE_based_gene/chr2_TE_based_gene_starts_count_win104")
gene_count=gcount[-1,3]
gene_count=as.vector(gene_count[-length(gene_count)])
pdf(file="MDS_1D_against_TE_based_gene_count_chr2_based_all_chr_win104.pdf")
plot(gene_count,x[(1+w[1]):sum(w[1:2])],ylab="Coordinate 1", main="MDS_1D_against_TE_based_gene_count_chr2_win104", pch=2, col=rainbow(8)[2])
dev.off()


gcount=read.table("~/Documents/Medicago_Project/gene_density_related/gene_and_TE_based_gene/chr3_TE_based_gene_starts_count_win104")
gene_count=gcount[-1,3]
gene_count=as.vector(gene_count[-length(gene_count)])
pdf(file="MDS_1D_against_TE_based_gene_count_chr3_based_all_chr_win104.pdf")
plot(gene_count,x[sum(w[1:2],1):sum(w[1:3])],ylab="Coordinate 1", main="MDS_1D_against_TE_based_gene_count_chr3_win104", pch=3, col=rainbow(8)[3])
dev.off()



gcount=read.table("~/Documents/Medicago_Project/gene_density_related/gene_and_TE_based_gene/chr4_TE_based_gene_starts_count_win104")
gene_count=gcount[-1,3]
gene_count=as.vector(gene_count[-length(gene_count)])
pdf(file="MDS_1D_against_TE_based_gene_count_chr4_based_all_chr_win104.pdf")
plot(gene_count,x[sum(w[1:3],1):sum(w[1:4])],ylab="Coordinate 1", main="MDS_1D_against_TE_based_gene_count_chr4_win104", pch=4, col=rainbow(8)[4])
dev.off()


gcount=read.table("~/Documents/Medicago_Project/gene_density_related/gene_and_TE_based_gene/chr5_TE_based_gene_starts_count_win104")
gene_count=gcount[-1,3]
gene_count=as.vector(gene_count[-length(gene_count)])
pdf(file="MDS_1D_against_TE_based_gene_count_chr5_based_all_chr_win104.pdf")
plot(gene_count,x[sum(w[1:4],1):sum(w[1:5])],ylab="Coordinate 1", main="MDS_1D_against_TE_based_gene_count_chr5_win104", pch=5, col=rainbow(8)[5])
dev.off()


gcount=read.table("~/Documents/Medicago_Project/gene_density_related/gene_and_TE_based_gene/chr6_TE_based_gene_starts_count_win104")
gene_count=gcount[-1,3]
gene_count=as.vector(gene_count[-length(gene_count)])
pdf(file="MDS_1D_against_TE_based_gene_count_chr6_based_all_chr_win104.pdf")
plot(gene_count,x[sum(w[1:5],1):sum(w[1:6])],ylab="Coordinate 1", main="MDS_1D_against_TE_based_gene_count_chr6_win104", pch=6, col=rainbow(8)[6])
dev.off()


gcount=read.table("~/Documents/Medicago_Project/gene_density_related/gene_and_TE_based_gene/chr7_TE_based_gene_starts_count_win104")
gene_count=gcount[-1,3]
gene_count=as.vector(gene_count[-length(gene_count)])
pdf(file="MDS_1D_against_TE_based_gene_count_chr7_based_all_chr_win104.pdf")
plot(gene_count,x[sum(w[1:6],1):sum(w[1:7])],ylab="Coordinate 1", main="MDS_1D_against_TE_based_gene_count_chr7_win104", pch=7, col=rainbow(8)[7])
dev.off()


gcount=read.table("~/Documents/Medicago_Project/gene_density_related/gene_and_TE_based_gene/chr8_TE_based_gene_starts_count_win104")
gene_count=gcount[-1,3]
gene_count=as.vector(gene_count[-length(gene_count)])
pdf(file="MDS_1D_against_TE_based_gene_count_chr8_based_all_chr_win104.pdf")
plot(gene_count,x[sum(w[1:7],1):sum(w[1:8])],ylab="Coordinate 1", main="MDS_1D_against_TE_based_gene_count_chr8_win104", pch=8, col=rainbow(8)[8])
dev.off()









