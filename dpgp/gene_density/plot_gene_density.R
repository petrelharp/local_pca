# downloaded 8/11/2016 from http://hgdownload.soe.ucsc.edu/goldenPath/dm3/database/refGene.txt.gz
# description from Table Scheme, RefGene, http://genome.ucsc.edu/cgi-bin/hgTables:
#
#   	Database: dm3    Primary Table: refGene    Row Count: 37,041   Data last updated: 2015-11-14
# Format description: A gene prediction with some additional info.
# field	example	SQL type 	info 	description
# bin 	614	smallint(5) unsigned 	range 	Indexing field to speed chromosome range queries.
# name 	NM_001316540	varchar(255) 	values 	Name of gene (usually transcript_id from GTF)
# chrom 	chr3R	varchar(255) 	values 	Reference sequence chromosome or scaffold
# strand 	-	char(1) 	values 	+ or - for strand
# txStart 	3856851	int(10) unsigned 	range 	Transcription start position
# txEnd 	3858731	int(10) unsigned 	range 	Transcription end position
# cdsStart 	3857933	int(10) unsigned 	range 	Coding region start
# cdsEnd 	3858629	int(10) unsigned 	range 	Coding region end
# exonCount 	1	int(10) unsigned 	range 	Number of exons
# exonStarts 	3856851,	longblob 	  	Exon start positions
# exonEnds 	3858731,	longblob 	  	Exon end positions
# score 	0	int(11) 	range 	score
# name2 	CG11035	varchar(255) 	values 	Alternate name (e.g. gene_id from GTF)
# cdsStartStat 	cmpl	enum('none', 'unk', 'incmpl', 'cmpl') 	values 	enum('none','unk','incmpl','cmpl')
# cdsEndStat 	cmpl	enum('none', 'unk', 'incmpl', 'cmpl') 	values 	enum('none','unk','incmpl','cmpl')
# exonFrames 	0,	longblob 	  	Exon frame {0,1,2}, or -1 if no frame for exon

lens <- structure(c(22963456L, 21142841L, 24536634L, 27894163L, 22418422L), .Dim = 5L, .Dimnames = list(c("chr2L", "chr2R", "chr3L", "chr3R", "chrX")))

genes <- read.table("refGene.txt.gz", header=FALSE, stringsAsFactors=FALSE)
names(genes) <- c("bin", "name", "chrom", "strand", "txStart", "txEnd", "cdsStart", 
            "cdsEnd", "exonCount", "exonStarts", "exonEnds", "score", "name2", 
            "cdsStartStat", "cdsEndStat", "exonFrames")
genes <- subset(genes, chrom %in% c("chr2L","chr2R","chr3L","chr3R","chrX"))
genes$chrom <- factor(genes$chrom)

window.length <- 1e5
gene.breaks <- seq(1,max(genes$txEnd)+window.length,by=window.length)
gene.mids <- gene.breaks[-1]-diff(gene.breaks)/2

# recomb rates
recomb <- read.table("recomb_rates.tsv",header=TRUE, sep='\t', stringsAsFactors=FALSE)
recomb$chrom <- paste0("chr",sapply( strsplit( recomb$Genomic.locus, ":" ), "[", 1 ))
recomb$start <- as.integer( sapply( strsplit( recomb$Genomic.locus, "[:.]" ), "[", 2 ) )
recomb$end <- as.integer( sapply( strsplit( recomb$Genomic.locus, "[:.]" ), "[", 4 ) )
recomb$mid <- (recomb$start+recomb$end)/2

# MDS scores: windows of 1,000 SNPs
mds.1e3 <- do.call( rbind, lapply( levels(genes$chrom), function (this.chrom) {
        read.table(file.path("..","results",sprintf("mds_%s_noinversion_window_1000snp.tsv",this.chrom)),sep='\t',header=TRUE)
    } ) )
mds.1e3$chrom <- gsub("Chr","chr",mds.1e3$chrom)
mds.1e3$mid <- (mds.1e3$start+mds.1e3$end)/2

# MDS scores: windows of 10,000 SNPs
mds.1e4 <- do.call( rbind, lapply( levels(genes$chrom), function (this.chrom) {
        read.table(file.path("..","results",sprintf("mds_%s_noinversion_window_10000snp.tsv",this.chrom)),sep='\t',header=TRUE)
    } ) )
mds.1e4$chrom <- gsub("Chr","chr",mds.1e4$chrom)
mds.1e4$mid <- (mds.1e4$start+mds.1e4$end)/2
mds.1e4$MDS1[ mds.1e4$chrom %in% c("chr3L","chr3R") ]  <- (-1) * mds.1e4$MDS1[ mds.1e4$chrom %in% c("chr3L","chr3R") ]

pdf( file="../../writeup/drosophila_recomb_mds.pdf", width=6, height=3, pointsize=10 )
layout( matrix( 1:(3*nlevels(genes$chrom)), nrow=3 ), 
       heights=c(1.25,1,1.5), widths=c(1.4,1,1,1,1) )
for (this.chrom in levels(genes$chrom)) {
    # png(file=sprintf("things_along_%s.png",this.chrom), width=4*144, height=6*144, res=144, pointsize=10)
    # layout(1:3,heights=c(1.2,1,1))
    m.left <- if ( this.chrom=="chr2L" ) { 4 } else { 0.25 }
    m.right <- if ( this.chrom=="chrX" ) { 1 } else { 0.25 }
    tx.bins <- as.numeric( table(cut(subset(genes,chrom==this.chrom)$txStart,breaks=gene.breaks)) 
                    + table(cut(subset(genes,chrom==this.chrom)$txEnd,breaks=gene.breaks)) )
    ylims <- c(0,1.1*quantile(tx.bins/2,.95))
    xlims <- c(0,lens[this.chrom]/1e6)
    par(mar=c(0.25,m.left,2,m.right)+.1)
    with( subset(mds.1e4,chrom==this.chrom), {
             plot( mid/1e6, MDS1, pch=20, col='red', cex=0.5,
                 xlim=xlims, main=this.chrom,
                 ylab=if(this.chrom=="chr2L"){'MDS 1'}else{""},
                 yaxt=if(this.chrom=="chr2L"){'s'}else{"n"},
                 xaxt='n', xlab='' )
        } )
    par(mar=c(0.25,m.left,0.25,m.right)+.1)
    with( subset(recomb,chrom==this.chrom), {
            plot( mid/1e6, (ylims[2]/15)*Comeron.Midpoint.rate, col='blue', pch=20,
                 xlim=xlims,cex=0.5,
               ylab=if(this.chrom=="chr2L"){'recombination rate'}else{""},
               yaxt=if(this.chrom=="chr2L"){'s'}else{"n"},
               xaxt='n', xlab='' );
        })
    par(mar=c(5,m.left,0.25,m.right)+.1)
    plot( gene.mids/1e6, tx.bins/2, cex=0.5,
           ylim=ylims, xlim=xlims, pch=20,
           ylab=if(this.chrom=="chr2L"){'gene count'}else{""},
           yaxt=if(this.chrom=="chr2L"){'s'}else{"n"},
           xlab='position (Mbp)' )
    # dev.off()
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

stats <- do.call( rbind, lapply( levels(genes$chrom), function (this.chrom) {
                this.start <- gene.breaks[which(gene.breaks < lens[this.chrom])]
                this.end <- gene.breaks[1+which(gene.breaks < lens[this.chrom])]
                this.breaks <- gene.breaks[c(1,1+which(gene.breaks < lens[this.chrom]))]
                data.frame(
                           chrom=this.chrom,
                           start=this.start,
                           end=this.end,
                           gene.count = with( subset(genes,chrom==this.chrom), ( as.numeric( table(cut(txStart,this.breaks)) ) + as.numeric( table(cut(txEnd,this.breaks)) ) )/2 ),
                           recomb = with( subset(recomb,chrom==this.chrom), average_windows( start, end, Comeron.Midpoint.rate, this.start, this.end ) ),
                           MDS1 = with( subset(mds.1e4,chrom==this.chrom), average_windows( start, end, MDS1, this.start, this.end ) )
                       )
           } ) )

tapply( 1:nrow(stats), stats$chrom, function (kk) {
           coef( summary( lm( recomb ~ MDS1, data=stats[kk,] ) ) )
        } )
# $chr2L
#             Estimate Std. Error   t value     Pr(>|t|)
# (Intercept) 2.633651  0.1221639 21.558336 4.684732e-56
# MDS1        5.916170  0.6602309  8.960759 1.436909e-16
# 
# $chr2R
#             Estimate Std. Error   t value     Pr(>|t|)
# (Intercept) 3.008740  0.1379687 21.807404 2.605543e-54
# MDS1        6.260817  0.9493917  6.594556 3.868908e-10
# 
# $chr3L
#              Estimate Std. Error   t value     Pr(>|t|)
# (Intercept)  1.977199 0.08985335 22.004731 1.736407e-59
# MDS1        -4.728210 0.59530596 -7.942488 7.648323e-14
# 
# $chr3R
#              Estimate Std. Error   t value     Pr(>|t|)
# (Intercept)  2.175621  0.1094316 19.881100 3.939434e-55
# MDS1        -4.867609  0.5725929 -8.500995 1.213363e-15
# 
# $chrX
#             Estimate Std. Error   t value     Pr(>|t|)
# (Intercept) 3.194880  0.1540236 20.742800 2.619987e-53
# MDS1        8.012976  0.9607161  8.340629 8.701074e-15

do.call( rbind, 
    tapply( 1:nrow(stats), stats$chrom, function (kk) {
               c( cor( stats$recomb[kk], stats$MDS1[kk], use='pairwise' ), 
               ( summary( lm( recomb ~ MDS1, data=stats[kk,] ) ) )$r.squared )
            } )
    )

# chr2L 0.5179585 0.2682810
# chr2R 0.4261314 0.1815880
# chr3L 0.4562214 0.2081380
# chr3R 0.4561819 0.2081019
# chrX  0.4935663 0.2436077

cor( stats$recomb, stats$MDS1, use='pairwise' )

pdf( file="../../writeup/drosophila_recomb_mds_correlation.pdf", width=6, height=2, pointsize=10 )
layout( t(1:5), widths=c(1.3,1,1,1,1.1) )
par( mar=c(5,4,2,0.25)+.1 )
for (this.chrom in levels(genes$chrom)) {
    with( subset( stats, chrom==this.chrom ), {
             plot( recomb, MDS1, pch=20, main=this.chrom,
                 xlab='recomb. rate',
                 yaxt=if(this.chrom=="chr2L") { 's' } else { 'n' },
                 ylab=if(this.chrom=="chr2L") { "MDS coordinate 1" } else { "" },
                  )
             abline( coef( lm( MDS1 ~ recomb ) ) )
        } )
    if (this.chrom=="chr3R") { par( mar=c(5,0.25,2,1)+.1 ) } else { par( mar=c(5,0.25,2,0.25)+.1 ) }
}
dev.off()



