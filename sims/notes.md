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
But this includes the centromere!  So, we're removing the windows overlapping with
anything on either side of the centromere out to the first segment in the genetic map having recombination rate at least 8cM/Mb,
which is from 13710990 to 16775980.

Note: the deCode map is obtained from http://www.decode.com/addendum/ ; should cite:
    Kong, A et al.  Fine scale recombination rate differences between sexes, populations and individuals. Nature  467 , 1099–1103 (28 October 2010) doi:10.1038/nature09525.
This file lists the position and the cM in the *previous* window, 
while msprime wants the cM/bp over the *next* window.


# Simulation of background selection:

For these we want to simulate a spatial population as above,
with flat recombination rate,
except in forwards time, 
and with 10,000 slightly deleterious loci 
(effect sizes drawn from an Exponential distribution with mean .001)
arranged with increasing density moving along the chromosome
(to mimic decreasing recombination rate).
Concretely, we let the inter-locus spacing be independent Exponentials, with the mean of the $k$th spacing
equal to $(1+9k/n)$, where $n$ is the number of loci;
and then we renormalize these positions to lie on the chromosome.
This implies density of selected loci will be ten times greater at the beginning of the chromosome than at the end.

Using the distribution of fitness effects from Eyre-Walker referenced by Kelley Harris: Gamma with shape=0.23 and mean=-0.043,
which translates to $\alpha=.23$ and $\beta=5.34$.
Output of .selloci is `loc a1 a2 fitness gen` (but appears not to have `gen`),
where `a1` and `a2` are alleles.

Modifying the above towards Drosophila:
all simulations on an 10x10 grid, with 2Ne=1000 per population,
a chromosome arm is about 25Mb,
mean recombination rate 2.5e-8, mutation rate 1e-7. 
Total population is 1e5, so we expect divergence of 0.01, and a tree every 400bp.
To have lineages mix over this time, but not instantly,
migration rate 4e-3 (in proportion of a population replaced by a given neighbor each generation).

Furthermore, we want to have some deleterious loci segregating:
a Drosophila chromosome arm is about 25Mb, of which about 1/5 is coding sequence;
divided into 1,000 loci this is 1Kb per locus;
if we assume that roughly 1/2 of these are deleterious
then we want a deleterious mutation rate of 5e-3.

Then, we want the universal ancestor to be maybe another 10,000 generations above the start of the simulation.

Timing results:
Here's the goal:
```
time ./background-sim.py -T 10000 -N 100 -w 10 -L 25e6 -l 1000 -m 4e-3 -u 5e-3 -r 2.5e-8 -a .23 -b 5.34 -s bground_sim_short.selloci \
            -A 10000 -k 1000 -U 1e-7 -o bground_sim_short.vcf -t bground_sim_short.trees -g bground_sim_short.log
```

```
# 32s : 100 gens, 2x2 grid, 10 samples, 100x smaller chromosome
time ./background-sim.py -T 100 -N 100 -w 2 -L .25e6 -l 1000 -m 4e-3 -u 5e-3 -r 2.5e-8 -a .23 -b 5.34 -s bground_sim_short.selloci \
            -A 10000 -k 10 -U 1e-7 -o bground_sim_short.vcf -t bground_sim_short.trees -g bground_sim_short.log

# 32s also : 100 gens, 2x2 grid, 100 samples, 100x smaller chromosome
time ./background-sim.py -T 100 -N 100 -w 2 -L .25e6 -l 1000 -m 4e-3 -u 5e-3 -r 2.5e-8 -a .23 -b 5.34 -s bground_sim_short.selloci \
            -A 10000 -k 100 -U 1e-7 -o bground_sim_short.vcf -t bground_sim_short.trees -g bground_sim_short.log

# 120s : 100 gens, 4x4 grid, 100 samples, 100x smaller chromosome
time ./background-sim.py -T 100 -N 100 -w 4 -L .25e6 -l 1000 -m 4e-3 -u 5e-3 -r 2.5e-8 -a .23 -b 5.34 -s bground_sim_short.selloci \
            -A 10000 -k 100 -U 1e-7 -o bground_sim_short.vcf -t bground_sim_short.trees -g bground_sim_short.log

# 150s: 100 gens, 4x4 grid, 1000 samples, 100x smaller chromosome
time ./background-sim.py -T 100 -N 100 -w 4 -L .25e6 -l 1000 -m 4e-3 -u 5e-3 -r 2.5e-8 -a .23 -b 5.34 -s bground_sim_short.selloci \
            -A 10000 -k 1000 -U 1e-7 -o bground_sim_short.vcf -t bground_sim_short.trees -g bground_sim_short.log

# 190 user/330 real: 100 gens, 4x4 grid, 1000 samples, 10x smaller chromosome
time ./background-sim.py -T 100 -N 100 -w 4 -L 2.5e6 -l 1000 -m 4e-3 -u 5e-3 -r 2.5e-8 -a .23 -b 5.34 -s bground_sim_short.selloci \
            -A 10000 -k 1000 -U 1e-7 -o bground_sim_short.vcf -t bground_sim_short.trees -g bground_sim_short.log

# 600s: 400 gens, 4x4 grid, 1000 samples, 10x smaller chromosome
# with vcf: 780/1260s but terminated writing 116G+ vcf file
time ./background-sim.py -T 400 -N 100 -w 4 -L 2.5e6 -l 1000 -m 4e-3 -u 5e-3 -r 2.5e-8 -a .23 -b 5.34 -s bground_sim_short.selloci \
            -A 10000 -k 1000 -U 1e-7 -o bground_sim_short.vcf -t bground_sim_short.trees -g bground_sim_short.log

# 40m sim + 10m simplify: 400 gens, 8x8 grid, 1000 samples
time ./background-sim.py -T 400 -N 100 -w 8 -L 25e6 -l 1000 -m 4e-3 -u 5e-3 -r 2.5e-8 -a .23 -b 5.34 -s bground_sim_short.selloci \
            -A 10000 -k 1000 -U 1e-7 -o bground_sim_short.vcf -t bground_sim_short.trees -g bground_sim_short.log

# 1000 gens, 8x8 grid
time ./background-sim.py -T 1000 -N 100 -w 8 -L 25e6 -l 1000 -m 4e-3 -u 5e-3 -r 2.5e-8 -a .23 -b 5.34 -s bground_sim_short.selloci \
            -A 10000 -k 1000 -U 1e-7 -o bground_sim_short.vcf -t bground_sim_short.trees -g bground_sim_short.log
```


With the following:
```
time ./background-sim.py -T 1000 -N 100 -w 8 -L 25e6 -l 1000 -m 4e-3 -u 5e-3 -r 2.5e-8 -a .23 -b 5.34 -s bground_sim_short.selloci \
            -A 10000 -k 1000 -U 1e-7 -o bground_sim_short.vcf -t bground_sim_short.trees -g bground_sim_short.log
```
Memory: 16Gb for 1000 gens.
2h to simulate; 1h to simplify.
