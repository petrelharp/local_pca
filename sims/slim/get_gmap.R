
## DeCode

gmap <- read.table('female.gmap',header=TRUE)
gmap <- subset(gmap, chr=="chr7" & !is.na(cM))
gmap <- subset(gmap, gmap$cM!=0 & c(gmap$cM[-1]==0,FALSE))
gmap$pos <- gmap$pos-gmap$pos[1]
gmap <- subset(gmap,c(TRUE,diff(gmap$pos)>0))
map <- data.frame(
    Chromosome=7,
    "Position(bp)"=gmap$pos,
    "Rate(cM/Mb)"=c(gmap$cM[-1]*1e6/diff(gmap$pos), 0.0),
    "Map(cM)"=0) # unused
# write.table(map, file="decode_female_chr7.gmap", row.names=FALSE)

# Chromosome  Position(bp)    Rate(cM/Mb)     Map(cM)
# chr1        55550   2.981822        0.000000
# chr1        82571   2.082414        0.080572

# SLiM wants endpoints and rates
write(gmap$pos[-1], file="decode_female_chr7.endpoints.txt")
write((gmap$cM[-1]/100)/diff(gmap$pos), file="decode_female_chr7.rates.txt")


## Step changes

relative_rates <- c(64, 16, 4, 1, 4, 16, 64)
step_endpoints <- floor(seq(0, max(gmap$pos), length.out=8))[-1]
step_rates <- sum(gmap$cM/100) * relative_rates / sum(relative_rates) / diff(c(0,step_endpoints))

write(step_endpoints, file="step_recomb.endpoints.txt")
write(step_rates, file="step_recomb.rates.txt")
