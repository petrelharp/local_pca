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

genes <- read.table("refGene.txt.gz", header=FALSE, stringsAsFactors=FALSE)
names(genes) <- c("bin", "name", "chrom", "strand", "txStart", "txEnd", "cdsStart", 
            "cdsEnd", "exonCount", "exonStarts", "exonEnds", "score", "name2", 
            "cdsStartStat", "cdsEndStat", "exonFrames")
genes <- subset(genes, chrom %in% c("chr2L","chr2R","chr3L","chr3R","chrX"))
genes$chrom <- factor(genes$chrom)

# recomb rates
recomb <- read.table("recomb_rates.tsv",header=TRUE, sep='\t', stringsAsFactors=FALSE)
recomb$chrom <- paste0("chr",sapply( strsplit( recomb$Genomic.locus, ":" ), "[", 1 ))
recomb$start <- as.integer( sapply( strsplit( recomb$Genomic.locus, "[:.]" ), "[", 2 ) )
recomb$end <- as.integer( sapply( strsplit( recomb$Genomic.locus, "[:.]" ), "[", 4 ) )
recomb$mid <- (recomb$start+recomb$end)/2

breaks <- seq(1,max(genes$txEnd)+1e5,by=5e5)
mids <- breaks[-1]-diff(breaks)/2

layout((1:nlevels(genes$chrom)))
par(mar=c(5,4,.5,4)+.1)
for (this.chrom in levels(genes$chrom)) {
    tx.bins <- as.numeric( table(cut(subset(genes,chrom==this.chrom)$txStart,breaks=breaks)) 
                    + table(cut(subset(genes,chrom==this.chrom)$txEnd,breaks=breaks)) )
    ylims <- c(0,1.1*quantile(tx.bins/2,.95))
    plot( mids,tx.bins/2,
            ylim=ylims, ylab=this.chrom, pch=20, cex=0.5 )
    with( subset(recomb,chrom==this.chrom), {
            points( mid, (ylims[2]/15)*Comeron.Midpoint.rate, col='red', pch=20 );
            axis(4, at=pretty(Comeron.Midpoint.rate)*(ylims[2]/15), 
                    labels=pretty(Comeron.Midpoint.rate), col='red', col.axis='red' )
            mtext("recombination rate", side=4, col='red', line=3, cex=0.75)
        })
}

