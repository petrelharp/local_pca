chr<-read.table(pipe("zcat POPRES_Genotypes_QC2_v2_TXT.tped.gz|grep '^1\\>'"))
geno <- as.matrix(chr[,-(1:4)])
geno[ ! (geno %in% c("A","C","G","T")) ] <- NA
Gs <- rowSums(geno=="G",na.rm=TRUE)
As <- rowSums(geno=="A",na.rm=TRUE)
Cs <- rowSums(geno=="C",na.rm=TRUE)
Ts <- rowSums(geno=="T",na.rm=TRUE)
triallelic <- ( ( (Gs==0) + (As==0) + (Cs==0) + (Ts==0) ) > 2 )
maxcounts <- pmax(Gs,As,Cs,Ts)
major <- rep(NA,nrow(geno))
major[ Gs==maxcounts ] <- "G"
major[ As==maxcounts ] <- "A"
major[ Ts==maxcounts ] <- "T"
major[ Cs==maxcounts ] <- "C"
stopifnot(sum(is.na(major))==0)
coded <- matrix( NA, nrow=nrow(geno), ncol(geno) )
coded[ geno==major ] <- 0
coded[ geno!=major ] <- 1
coded[is.na(geno)]<-NA
n<-ncol(coded)/2
coded<-coded[,2*(1:n)]+coded[,2*(1:n)-1]
write.table(coded,"coded_data_chr1.txt",sep="\t")
