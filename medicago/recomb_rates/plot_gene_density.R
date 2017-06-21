# using data from Tim Paape, 8/17/2016:
# Hi Peter R
# 
# Here are the complete LDhat output files for all 1kb sliding window and all
# chromosomes, with positions of each window on each chromosome (sliding
# windows).  In our GBE paper we published only chromosomes 2, 3 and 5. The
# other 5 chromosomes are unpublished recombination rates  for Medicago, so
# depending on what you want to do perhaps we can have that discussion later. 
# 
# This data comes from SNPs called using our version 3.0 gene annotations (also
# used in the Branca et al. 2011 PNAS paper, 26 accessions). Since then there
# have been v3.5 and 4.0 gene annotations and possible slightly modified
# scaffolds in the genome assembly versions.  Also, if you need I can find the
# centromere positions if you need them for v3.0.
# 
# cheers
# Tim

recomb <- read.table("from_paape/07_recom.txt",header=TRUE)
recomb$mid <- (recomb$sStart + recomb$sEnd)/2
levels(recomb$chr) <- gsub("MtChr","chr",levels(recomb$chr))

# map <- read.csv("from_paape/Mtruncatula_mapbased_recomb_2may_b.csv",header=TRUE)
# map <- map[,c("chr","start","end","cM","cm.bp...1.000.000","ave.cM...Mbp....3.window.moving.average")]

mds <- do.call( rbind, lapply( 1:8, function (k) {
                          read.table(file.path("..","results",sprintf("mds_chr%d_window_10000snp.tsv",k)), header=TRUE)
        } ) )
mds$mid <- (mds$start+mds$end)/2
levels(mds$chrom) <- tolower(levels(mds$chrom))

xlims <- lapply(levels(mds$chrom), function (this.chrom) { range(subset(recomb,chr==this.chrom)$mid, subset(mds,chrom==this.chrom)$mid) } )
names(xlims) <- levels(mds$chrom)

pdf(file="recomb_and_mds.pdf", width=6, height=6, pointsize=10)
layout( matrix(1:16,nrow=4), heights=c(1,1.2,1,1.2), widths=c(1.2,1,1,1) )
for (this.chrom in paste0('chr',1:8)) {
    left.mar <- if (this.chrom %in% paste0('chr',1:2)) { 4} else {0}
    par(mar=c(0.25,left.mar,2,2)+.1)
    with( subset(recomb,chr==this.chrom), {
         plot(mid, rho.mean./1000, pch=20,
              xlim=xlims[[this.chrom]],
              xaxt='n',
              ylab='recomb rate (cM/Mb)')
         mtext(3,text=this.chrom)
    } )
    par(mar=c(5,left.mar,0.25,2)+.1)
    with( subset(mds,chrom==this.chrom),
         plot(mid,MDS1,col='red',
              xlim=xlims[[this.chrom]],
              ylab='MDS 1',
              xlab='position (bp)',
              pch=20)
     )
}
dev.off()

average_windows <- function (start,end,value,new.start,new.end) {
    out <- numeric(length(new.start))
    for (k in seq_along(new.start)) {
        weights <- pmax(0,pmin( new.end[k], end ) - pmax( new.start[k], start )) / (end-start)
        out[k] <- sum(weights * value)/sum(weights)
    }
    return(out)
}

stats <- do.call( rbind, lapply( levels(mds$chrom), function (this.chrom) {
                this.breaks <- with(subset(mds,chrom==this.chrom), 
                                    c( start[1], (start[-1]+end[-length(end)])/2, end[length(end)] ) 
                                )
                this.start <- this.breaks[-length(this.breaks)]
                this.end <- this.breaks[-1]
                data.frame(
                           chrom=this.chrom,
                           start=this.start,
                           end=this.end,
                           recomb = with( subset(recomb,chr==this.chrom), 
                                         average_windows( sStart, sEnd, rho.mean., this.start, this.end ) ),
                           MDS1 = with( subset(mds,chrom==this.chrom), 
                                       average_windows( start, end, MDS1, this.start, this.end ) )
                       )
           } ) )

tapply( 1:nrow(stats), stats$chrom, function (kk) {
           coef( summary( lm( recomb ~ MDS1, data=stats[kk,] ) ) )
        } )

do.call( rbind, 
    tapply( 1:nrow(stats), stats$chrom, function (kk) {
               c( cor=cor( stats$recomb[kk], stats$MDS1[kk], use='pairwise' ), 
               r.squared=( summary( lm( recomb ~ MDS1, data=stats[kk,] ) ) )$r.squared )
            } )
    )

layout( matrix(1:8, nrow=4), widths=c(1.3,1,1,1,1.1) )
par( mar=c(5,4,2,0.25)+.1 )
for (this.chrom in levels(mds$chrom)) {
    with( subset( stats, chrom==this.chrom ), {
             plot( recomb, MDS1, pch=20, main=this.chrom,
                 xlab='recomb. rate',
                 yaxt=if(this.chrom=="chr2L") { 's' } else { 'n' },
                 ylab=if(this.chrom=="chr2L") { "MDS coordinate 1" } else { "" },
                  )
             abline( coef( lm( MDS1 ~ recomb ) ) )
        } )
    if (this.chrom=="chr7") { par( mar=c(5,0.25,2,1)+.1 ) } else { par( mar=c(5,0.25,2,0.25)+.1 ) }
}


