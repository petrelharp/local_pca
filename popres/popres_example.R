library(lostruct)

# as in POPRES_SNPdata_recode12.R
chr22 <- read_tped("POPRES_Genotypes_QC2_v2_TXT.tped.gz", chrom=22)

# as in POPRES_PCA_win100.R
eigenstuff <- eigen_windows(chr22, win=100, k=2)

# as in POPRES_distance.R
windist <- pc_dist( eigenstuff, npc=2 )

# as in POPRES_MDS.R
fit2d <- cmdscale( windist, eig=TRUE, k=2 )
plot( fit2d$points, xlab="Coordinate 1", ylab="Coordinate 2", col=rainbow(1.2*nrow(windist)) )
