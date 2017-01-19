# Neutral simulations:

In `neutral/`; all simulations on an 10x10 grid, with 2Ne=1000 per population,
simulating a chromosome of length 4e7bp, mean recombination rate 2.5e-8, mutation rate 1e-7. 
Total population is 1e5, so we expect divergence of 0.01, and a tree every 400bp.
To have lineages mix over this time, but not instantly,
migration rate 4e-3 (in proportion of a population replaced by a given neighbor each generation).

In all cases we sample 1000 individuals.

1. Flat recombination:
```
./neutral-sim.py -k 36   -N 100  -w 3  -L 4e7 -m 4e-3 -r 2.5e-8 -R 2.5e-8 -u 1e-7 -o flat_recomb_short.vcf -g flat_recomb_short.log
./neutral-sim.py -k 1000 -N 1000 -w 10 -L 4e7 -m 4e-3 -r 2.5e-8 -R 2.5e-8 -u 1e-7 -o flat_recomb.vcf -g flat_recomb.log
bcftools convert -O b -o flat_recomb.bcf flat_recomb.vcf
```

2. Linearly increasing recombination rate from 0.5e-8 to 4.5e-8 (same total length):
```
./neutral-sim.py -k 36   -N 100  -w 3  -L 4e7 -m 4e-3 -r 0.5e-8 -R 4.5e-8 -u 1e-7 -o linear_recomb_short.vcf -g linear_recomb_short.log
./neutral-sim.py -k 1000 -N 1000 -w 10 -L 4e7 -m 4e-3 -r 0.5e-8 -R 4.5e-8 -u 1e-7 -o linear_recomb.vcf -g linear_recomb.log
bcftools convert -O b -o linear_recomb.bcf linear_recomb.vcf
```

3. Hot spottish recombination map (same total length):
The map we use is deCode's female map for chr7, rescaled to our length:
```r
gmap <- read.table('female.gmap',header=TRUE)
gmap <- subset(gmap, chr=="chr7" & !is.na(cM))
gmap <- subset(gmap, gmap$cM!=0 & c(gmap$cM[-1]==0,FALSE))
bp_len <- 4e7
gmap$pos <- gmap$pos-gmap$pos[1]
gmap$pos <- floor( gmap$pos*bp_len/max(gmap$pos) )
gmap$cM <- gmap$cM * 100 / sum(gmap$cM)
gmap <- subset(gmap,c(TRUE,diff(gmap$pos)>0))
map <- data.frame(
    Chromosome=1,
    "Position(bp)"=gmap$pos,
    "Rate(cM/Mb)"=c( gmap$cM[-1]*1e6/diff(gmap$pos), 0.0),
    "Map(cM)"=0) # unused
write.table(map, file="decode_female_chr7.gmap", row.names=FALSE)

# Chromosome  Position(bp)    Rate(cM/Mb)     Map(cM)
# chr1        55550   2.981822        0.000000
# chr1        82571   2.082414        0.080572
```
```
./neutral-sim.py -k 36   -N 100  -w 3  -L 4e7 -m 4e-3 -p decode_female_chr7.gmap -u 1e-7 -o hot_recomb_short.vcf -g hot_recomb_short.log
./neutral-sim.py -k 1000 -N 1000 -w 10 -L 4e7 -m 4e-3 -p decode_female_chr7.gmap -u 1e-7 -o hot_recomb.vcf -g hot_recomb.log
bcftools convert -O b -o hot_recomb.bcf hot_recomb.vcf
```

Note: the deCode map is obtained from http://www.decode.com/addendum/ ; should cite:
    Kong, A et al.  Fine scale recombination rate differences between sexes, populations and individuals. Nature  467 , 1099–1103 (28 October 2010) doi:10.1038/nature09525.
This file lists the position and the cM in the *previous* window, 
while msprime wants the cM/bp over the *next* window.


# Simulation of background selection:

For these we want to simulate a spatial population as above,
with flat recombination rate,
except in forwrds time, 
and with 10,000 slightly deleterious loci 
(effect sizes drawn from an Exponential distribution with mean .001)
arranged with increasing density moving along the chromosome
(to mimic decreasing recombination rate).
Concretely, we let the inter-locus spacing be independent Exponentials, with the mean of the $k$th spacing
equal to $(1+9k/n)$, where $n$ is the number of loci;
and then we renormalize these positions to lie on the chromosome.
This implies density of selected loci will be ten times greater at the beginning of the chromosome than at the end.

