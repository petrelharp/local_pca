#/bin/env Rscript

usage <- "
Gets some summary statistics comparing the results (MDS coordinates)
across runs of the algorithm with different parameters.
"

library(jsonlite)

dirs <- c("./lostruct_results/type_snp_size_10000_weights_none_jobid_278544", 
            "./lostruct_results/type_snp_size_1000_weights_none_jobid_450751", 
            "./lostruct_results/type_snp_size_10000_weights_none_jobid_080290", 
            "./lostruct_results/type_bp_size_100000_weights_none_jobid_381845", 
            "./lostruct_results/type_bp_size_10000_weights_none_jobid_519007")

param.list <- lapply(dirs, function (dd) {
            fromJSON(file.path(dd,"config.json")) } )

params <- data.frame(
        outdir=sapply(param.list,"[[","outdir"),
        type=sapply(param.list,"[[","type"),
        size=sapply(param.list,"[[","size"),
        npc=sapply(param.list,"[[","npc"),
        nmds=sapply(param.list,"[[","nmds")
        )

chrom.names <- paste0("chr",1:8)
region.file.names <- paste0(chrom.names,"-filtered-set-2014Apr15.regions.csv")
chrom.lens <- c( chr1=52991155, chr2=45729672, chr3=55515152, chr4=56582383, chr5=43630510, chr6=35275713, chr7=49172423, chr8=45569985, chl_Mt=124033 )
chrom.starts <- cumsum(c(0,chrom.lens[-length(chrom.lens)]))
names(chrom.starts) <- names(chrom.lens)
chrom_pos <- function (chrom,pos) {
    return(pos + chrom.starts[chrom])
}

get_regions <- function (dd) {
    out <- do.call(rbind, lapply( file.path(dd,region.file.names), function (dn) {
                        z <- read.csv(dn, header=TRUE, stringsAsFactors=FALSE)
                        this.chrom <- z$chrom[1]
                        breaks <- c(0,(1/2)*(z$start[-1]+z$end[-nrow(z)]),chrom.lens[this.chrom])
                        z$real_start <- chrom.starts[this.chrom]+breaks[-nrow(z)]
                        z$real_end <- chrom.starts[this.chrom]+breaks[-1]
                        return(z)
                    } ) )
}

match_window <- function (chrom,pos,reg) {
    # find which window corresp to (chrom,pos) in reg
    cp <- chrom_pos(chrom,pos)
    return(findInterval(cp,c(0,reg$real_end)))
}

compare_mds <- function (d1,d2,k) {
    # correlation of d1 with mean of matching windows in d2
    reg1 <- get_regions(d1)
    reg2 <- get_regions(d2)
    win2 <- factor(match_window(reg2$chrom,(reg2$start+reg2$end)/2,reg1),levels=1:nrow(reg1))
    mds1 <- read.csv(file.path(d1,"mds_coords.csv"),header=TRUE)
    mds2 <- read.csv(file.path(d2,"mds_coords.csv"),header=TRUE)
    nmds1 <- sum(grepl("MDS",colnames(mds1)))
    nmds2 <- sum(grepl("MDS",colnames(mds2)))
    nmds <- min(nmds1,nmds2)
    if (k>nmds) { return(NA) }
    this.mds2 <- tapply(mds2[,paste0("MDS",k)],win2,mean,na.rm=TRUE)
    return( cor( mds1[,paste0("MDS",k)], this.mds2, use="pairwise" ) )
           
}

# Produces a matrix with upper triangle correlations in MDS1, and lower triangle in MDS2
mds.cors <- list( matrix(NA,nrow=nrow(params),ncol=ncol(params)) )[c(1,1)]
for (i in 1:nrow(params)) {
    for (j in 1:nrow(params)) {
        mds.cors[[1]][i,j] <- compare_mds(params$outdir[i],params$outdir[j],1)
        mds.cors[[2]][j,i] <- compare_mds(params$outdir[i],params$outdir[j],2)
    }
}
for (k in 1:2) {
    colnames(mds.cors[[k]]) <- rownames(mds.cors[[k]]) <- sprintf("%d%s, %d PCs", params$size, params$type, params$npc)
}

library(xtable)

options(digits=2)
lapply(mds.cors,xtable)

# % latex table generated in R 3.3.1 by xtable 1.8-2 package
# % Wed Feb  8 16:17:12 2017
# \begin{table}[ht]
# \centering
# \begin{tabular}{rrrrrr}
#   \hline
#  & 10000snp, 2 PCs & 1000snp, 2 PCs & 10000snp, 5 PCs & 100000bp, 2 PCs & 10000bp, 2 PCs \\ 
#   \hline
# 10000snp, 2 PCs & 1.00 & 0.87 & 0.96 & 0.90 & 0.88 \\ 
#   1000snp, 2 PCs & 0.68 & 1.00 & 0.73 & 0.68 & 0.94 \\ 
#   10000snp, 5 PCs & 0.96 & 0.92 & 1.00 & 0.88 & 0.93 \\ 
#   100000bp, 2 PCs & 0.90 & 0.87 & 0.88 & 1.00 & 0.87 \\ 
#   10000bp, 2 PCs & 0.68 & 0.93 & 0.72 & 0.67 & 1.00 \\ 
#    \hline
# \end{tabular}
# \end{table}
# 
# [[2]]
# % latex table generated in R 3.3.1 by xtable 1.8-2 package
# % Wed Feb  8 16:17:12 2017
# \begin{table}[ht]
# \centering
# \begin{tabular}{rrrrrr}
#   \hline
#  & 10000snp, 2 PCs & 1000snp, 2 PCs & 10000snp, 5 PCs & 100000bp, 2 PCs & 10000bp, 2 PCs \\ 
#   \hline
# 10000snp, 2 PCs & 1.00 & 0.54 & 0.93 & 0.87 & 0.56 \\ 
#   1000snp, 2 PCs & 0.82 & 1.00 & 0.76 & 0.83 & 0.92 \\ 
#   10000snp, 5 PCs & 0.93 & 0.50 & 1.00 & 0.83 & 0.52 \\ 
#   100000bp, 2 PCs & 0.87 & 0.59 & 0.84 & 1.00 & 0.58 \\ 
#   10000bp, 2 PCs & 0.83 & 0.92 & 0.77 & 0.84 & 1.00 \\ 
#    \hline
# \end{tabular}
# \end{table}
# 
