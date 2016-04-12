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

gcount.files <- c(
        "~/Documents/Medicago_Project/gene_density_related/gene_and_TE_based_gene/chr1_gene_starts_count_win104",
        "~/Documents/Medicago_Project/gene_density_related/gene_and_TE_based_gene/chr2_gene_starts_count_win104",
        "~/Documents/Medicago_Project/gene_density_related/gene_and_TE_based_gene/chr3_gene_starts_count_win104",
        "~/Documents/Medicago_Project/gene_density_related/gene_and_TE_based_gene/chr4_gene_starts_count_win104",
        "~/Documents/Medicago_Project/gene_density_related/gene_and_TE_based_gene/chr5_gene_starts_count_win104",
        "~/Documents/Medicago_Project/gene_density_related/gene_and_TE_based_gene/chr6_gene_starts_count_win104",
        "~/Documents/Medicago_Project/gene_density_related/gene_and_TE_based_gene/chr7_gene_starts_count_win104",
        "~/Documents/Medicago_Project/gene_density_related/gene_and_TE_based_gene/chr8_gene_starts_count_win104"
    )

plot.gcount <- function (chrom) {
    gcount=read.table(gcount.files[chrom])
    gene_count=gcount[-1,3]
    gene_count=as.numeric(as.vector(gene_count[-length(gene_count)]))
    sapply(1:w[chrom], function(i){(snp_pos[i*10000]+snp_pos[(i-1)*10000+1])/2}) -> mid_pos
    floor(mid_pos) -> wid_mid_pos
    outfile <- sprintf("MDS_1D_and_gene_count_and_lm_chr%d_based_all_chr_win104.pdf",chrom)
    pdf(file="MDS_1D_and_gene_count_and_lm_chr1_based_all_chr_win104.pdf",width=6,height=3, pointsize=10)
    layout(matrix(c(1,2,3,3), nrow=2),heights=c(2,2,4))
    par(mar=c(4.5,5,1.5,1))
    plot(wid_mid_pos,x[1:w[1]],ylab="Coordinate 1", main=sprintf("MDS vs position, chr. %d", chrom), pch=1, col=rainbow(8)[1])
    plot(wid_mid_pos,gene_count)
    plot(gene_count,x[1:w[1]],ylab="Coordinate 1", main=sprintf("MDS vs gene count, chr. %d", chrom), pch=1, col=rainbow(8)[1])
    regl<-lm(x[1:w[chrom]]~gene_count)
    abline(regl)
    dev.off()
}

for (k in seq_along(gcount.files)) { plot.gcount(k) }

# 
# gcount=read.table("~/Documents/Medicago_Project/gene_density_related/gene_and_TE_based_gene/chr1_gene_starts_count_win104")
# gene_count=gcount[-1,3]
# gene_count=as.numeric(as.vector(gene_count[-length(gene_count)]))
# sapply(1:n1, function(i){(snp_pos[i*10000]+snp_pos[(i-1)*10000+1])/2}) -> mid_pos
# floor(mid_pos) -> wid_mid_pos
# pdf(file="MDS_1D_and_gene_count_and_lm_chr1_based_all_chr_win104.pdf",width=6,height=3, pointsize=10)
# layout(matrix(c(1,2,3,3), nrow=2),heights=c(2,2,4))
# par(mar=c(4.5,5,1.5,1))
# plot(wid_mid_pos,x[1:w[1]],ylab="Coordinate 1", main="MDS_1D_and_gene_count_chr1_win104", pch=1, col=rainbow(8)[1])
# plot(wid_mid_pos,gene_count)
# plot(gene_count,x[1:w[1]],ylab="Coordinate 1", main="MDS_1D_against_gene_count_chr1_win104_with_lm", pch=1, col=rainbow(8)[1])
# regl<-lm(x[1:w[1]]~gene_count)
# abline(regl)
# dev.off()
# 
# 
# gcount=read.table("~/Documents/Medicago_Project/gene_density_related/gene_and_TE_based_gene/chr2_gene_starts_count_win104")
# gene_count=gcount[-1,3]
# gene_count=as.numeric(as.vector(gene_count[-length(gene_count)]))
# sapply(1:n2, function(i){(snp_pos[i*10000]+snp_pos[(i-1)*10000+1])/2}) -> mid_pos
# floor(mid_pos) -> wid_mid_pos
# pdf(file="MDS_1D_and_gene_count_and_lm_chr2_based_all_chr_win104.pdf",width=6,height=3,pointsize=10)
# layout(matrix(c(1,2,3,3), nrow=2),heights=c(2,2,4))
# par(mar=c(4.5,5,1.5,1))
# plot(wid_mid_pos,x[(1+w[1]):sum(w[1:2])],ylab="Coordinate 1", main="MDS_1D_and_gene_count_chr2_win104", pch=1, col=rainbow(8)[1])
# plot(wid_mid_pos,gene_count)
# plot(gene_count,x[(1+w[1]):sum(w[1:2])],ylab="Coordinate 1", main="MDS_1D_against_gene_count_chr2_win104_with_lm", pch=2, col=rainbow(8)[2])
# regl<-lm(x[(1+w[1]):sum(w[1:2])]~gene_count)
# abline(regl)
# dev.off()
# 
# gcount=read.table("~/Documents/Medicago_Project/gene_density_related/gene_and_TE_based_gene/chr3_gene_starts_count_win104")
# gene_count=gcount[-1,3]
# gene_count=as.numeric(as.vector(gene_count[-length(gene_count)]))
# sapply(1:n3, function(i){(snp_pos[i*10000]+snp_pos[(i-1)*10000+1])/2}) -> mid_pos
# floor(mid_pos) -> wid_mid_pos
# pdf(file="MDS_1D_and_gene_count_and_lm_chr3_based_all_chr_win104.pdf",width=6,height=3)
# layout(matrix(c(1,2,3,3), nrow=2),heights=c(2,2,4))
# par(mar=c(4.5,5,1.5,1))
# plot(wid_mid_pos,x[sum(w[1:2],1):sum(w[1:3])],ylab="Coordinate 1", main="MDS_1D_and_gene_count_chr3_win104", pch=1, col=rainbow(8)[1])
# plot(wid_mid_pos,gene_count)
# plot(gene_count,x[sum(w[1:2],1):sum(w[1:3])],ylab="Coordinate 1", main="MDS_1D_against_gene_count_chr3_win104_with_lm", pch=3, col=rainbow(8)[3])
# regl<-lm(x[sum(w[1:2],1):sum(w[1:3])]~gene_count)
# abline(regl)
# dev.off()
# 
# 
# gcount=read.table("~/Documents/Medicago_Project/gene_density_related/gene_and_TE_based_gene/chr4_gene_starts_count_win104")
# gene_count=gcount[-1,3]
# gene_count=as.numeric(as.vector(gene_count[-length(gene_count)]))
# sapply(1:n4, function(i){(snp_pos[i*10000]+snp_pos[(i-1)*10000+1])/2}) -> mid_pos
# floor(mid_pos) -> wid_mid_pos
# pdf(file="MDS_1D_and_gene_count_and_lm_chr4_based_all_chr_win104.pdf",width=6,height=3)
# layout(matrix(c(1,2,3,3), nrow=2),heights=c(2,2,4))
# par(mar=c(4.5,5,1.5,1))
# plot(wid_mid_pos,x[sum(w[1:3],1):sum(w[1:4])],ylab="Coordinate 1", main="MDS_1D_and_gene_count_chr4_win104", pch=1, col=rainbow(8)[1])
# plot(wid_mid_pos,gene_count)
# plot(gene_count,x[sum(w[1:3],1):sum(w[1:4])],ylab="Coordinate 1", main="MDS_1D_against_gene_count_chr4_win104_with_lm", pch=4, col=rainbow(8)[4])
# regl<-lm(x[sum(w[1:3],1):sum(w[1:4])]~gene_count)
# abline(regl)
# dev.off()
# 
# gcount=read.table("~/Documents/Medicago_Project/gene_density_related/gene_and_TE_based_gene/chr5_gene_starts_count_win104")
# gene_count=gcount[-1,3]
# gene_count=as.numeric(as.vector(gene_count[-length(gene_count)]))
# sapply(1:n5, function(i){(snp_pos[i*10000]+snp_pos[(i-1)*10000+1])/2}) -> mid_pos
# floor(mid_pos) -> wid_mid_pos
# pdf(file="MDS_1D_and_gene_count_and_lm_chr5_based_all_chr_win104.pdf",width=10,height=5)
# layout(matrix(c(1,2,3,3), nrow=2),heights=c(2,2,4))
# par(mar=c(4.5,5,1.5,1))
# plot(wid_mid_pos,x[sum(w[1:4],1):sum(w[1:5])],ylab="Coordinate 1", main="MDS_1D_and_gene_count_chr5_win104", pch=1, col=rainbow(8)[1])
# plot(wid_mid_pos,gene_count)
# plot(gene_count,x[sum(w[1:4],1):sum(w[1:5])],ylab="Coordinate 1", main="MDS_1D_against_gene_count_chr5_win104_with_lm", pch=5, col=rainbow(8)[5])
# regl<-lm(x[sum(w[1:4],1):sum(w[1:5])]~gene_count)
# abline(regl)
# dev.off()
# 
# gcount=read.table("~/Documents/Medicago_Project/gene_density_related/gene_and_TE_based_gene/chr6_gene_starts_count_win104")
# gene_count=gcount[-1,3]
# gene_count=as.numeric(as.vector(gene_count[-length(gene_count)]))
# sapply(1:n6, function(i){(snp_pos[i*10000]+snp_pos[(i-1)*10000+1])/2}) -> mid_pos
# floor(mid_pos) -> wid_mid_pos
# pdf(file="MDS_1D_and_gene_count_and_lm_chr6_based_all_chr_win104.pdf",width=10,height=5)
# layout(matrix(c(1,2,3,3), nrow=2),heights=c(2,2,4))
# par(mar=c(4.5,5,1.5,1))
# plot(wid_mid_pos,x[sum(w[1:5],1):sum(w[1:6])],ylab="Coordinate 1", main="MDS_1D_and_TE_based_gene_count_chr6_win104", pch=1, col=rainbow(8)[1])
# plot(wid_mid_pos,gene_count)
# plot(gene_count,x[sum(w[1:5],1):sum(w[1:6])],ylab="Coordinate 1", main="MDS_1D_against_gene_count_chr6_win104_with_lm", pch=6, col=rainbow(8)[6])
# regl<-lm(x[sum(w[1:5],1):sum(w[1:6])]~gene_count)
# abline(regl)
# dev.off()
# 
# 
# gcount=read.table("~/Documents/Medicago_Project/gene_density_related/gene_and_TE_based_gene/chr7_gene_starts_count_win104")
# gene_count=gcount[-1,3]
# gene_count=as.numeric(as.vector(gene_count[-length(gene_count)]))
# sapply(1:n7, function(i){(snp_pos[i*10000]+snp_pos[(i-1)*10000+1])/2}) -> mid_pos
# floor(mid_pos) -> wid_mid_pos
# pdf(file="MDS_1D_and_gene_count_and_lm_chr7_based_all_chr_win104.pdf",width=10,height=5)
# layout(matrix(c(1,2,3,3), nrow=2),heights=c(2,2,4))
# par(mar=c(4.5,5,1.5,1))
# plot(wid_mid_pos,x[sum(w[1:6],1):sum(w[1:7])],ylab="Coordinate 1", main="MDS_1D_and_gene_count_chr7_win104", pch=1, col=rainbow(8)[1])
# plot(wid_mid_pos,gene_count)
# plot(gene_count,x[sum(w[1:6],1):sum(w[1:7])],ylab="Coordinate 1", main="MDS_1D_against_gene_count_chr7_win104_with_lm", pch=7, col=rainbow(8)[7])
# regl<-lm(x[sum(w[1:6],1):sum(w[1:7])]~gene_count)
# abline(regl)
# dev.off()
# 
# 
# gcount=read.table("~/Documents/Medicago_Project/gene_density_related/gene_and_TE_based_gene/chr8_gene_starts_count_win104")
# gene_count=gcount[-1,3]
# gene_count=as.numeric(as.vector(gene_count[-length(gene_count)]))
# sapply(1:n8, function(i){(snp_pos[i*10000]+snp_pos[(i-1)*10000+1])/2}) -> mid_pos
# floor(mid_pos) -> wid_mid_pos
# pdf(file="MDS_1D_and_gene_count_and_lm_chr8_based_all_chr_win104.pdf",width=10,height=5)
# layout(matrix(c(1,2,3,3), nrow=2),heights=c(2,2,4))
# par(mar=c(4.5,5,1.5,1))
# plot(wid_mid_pos,x[sum(w[1:7],1):sum(w[1:8])],ylab="Coordinate 1", main="MDS_1D_and_gene_count_chr8_win104", pch=1, col=rainbow(8)[1])
# plot(wid_mid_pos,gene_count)
# plot(gene_count,x[sum(w[1:7],1):sum(w[1:8])],ylab="Coordinate 1", main="MDS_1D_against_gene_count_chr8_win104_with_lm", pch=8, col=rainbow(8)[8])
# regl<-lm(x[sum(w[1:7],1):sum(w[1:8])]~gene_count)
# abline(regl)
# dev.off()
