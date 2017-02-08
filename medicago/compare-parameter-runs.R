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
mds.cors <- matrix(NA,nrow=nrow(params),ncol=ncol(params))
for (i in 1:(nrow(params)-1)) {
    for (j in (i+1):nrow(params)) {
        mds.cors[i,j] <- compare_mds(params$outdir[i],params$outdir[j],1)
        mds.cors[j,i] <- compare_mds(params$outdir[i],params$outdir[j],2)
    }
}
rownames(mds.cors) <- sprintf("window=%d%s, %d PCs", params$size, params$type, params$npc)
mds.cors

#                             [,1]      [,2]      [,3]      [,4]      [,5]
# window=10000snp, 2 PCs        NA 0.8701600 0.9593167 0.9024207 0.8765862
# window=1000snp, 2 PCs  0.8195672        NA 0.7253759 0.6820970 0.9355702
# window=10000snp, 5 PCs 0.9319838 0.5027238        NA 0.8831645 0.9284340
# window=100000bp, 2 PCs 0.8724878 0.5896812 0.8378923        NA 0.8724009
# window=10000bp, 2 PCs  0.8279086 0.9172536 0.7738415 0.8373688        NA

