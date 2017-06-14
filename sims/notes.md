# Recap of the below

- [background, 5x5 grid](background/bground_sim_5000gens_5x5_migr0.01_recomb_1e-7/type_snp_size_200_npc_2_jobid_001/run-summary.html) -
    intriguing noisy patterns that don't line up with variation in recombination rate
- [background, 5x5 grid](background/bground_sim_5000gens_5x5/type_snp_size_500_npc_2_jobid_003/run_summary.html) -
    with a different migraiton rate to the above; intriguing smooth patterns in spikes along the chromosome
- [symmetric local](local/threeway_sym_200_1000_0.005_10.0_.000004_10_19024/bp_20000_npc_2/run-summary.html) -
    three pops with alternate halves of the chromosome carrying differently locally adaptive alleles
    three corners reflect greater or lesser spreading out of populations,
    somewhat partitioned in halves but noisy
    * [more windows](local/threeway_sym_200_1000_0.005_10.0_.000004_10_19024/mu_1e-4_bp_2000_npc_2/run-summary.html) same thing with higher mutation rate and more windows

**Pseudo-simulations:**

- [threesim](background/threesim_1e6_8657/bp_10000_npc_2/run-summary.html) - 
    separate, neutral simulations of two chromosomes mimicking a
    recent switch of `A|B->C` to `A<-B|C` - has good separation on MDS 1
- [threesim](background/threesim_1e6_16001/bp_10000_npc_2/run-summary.html) -
    another example


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

Modifying the above towards Drosophila: if possible,
all simulations on an 10x10 grid, with 2Ne=1000 per population,
a chromosome arm is about 25Mb,
neutral mutation rate 1e-7. 
Recombination rate in Drosophila is around 2.5 cM/Mb, or 2.5e-8 per bp,
which will give 0.625 recombinations per chromosome per generation.
Total population is 1e5, so we expect divergence of 0.01, and a tree every 400bp.
To have lineages mix over this time, but not instantly,
migration rate 4e-3 (in proportion of a population replaced by a given neighbor each generation).

Furthermore, we want to have some deleterious loci segregating:
a Drosophila chromosome arm is about 25Mb, of which about 1/5 is coding sequence;
divided into 1,000 loci this is 1Kb per locus;
if we assume that roughly 1/2 of these are deleterious
then we want a deleterious mutation rate of 5e-3.
This has a total mutation rate of 5 per generation, or 0.2e-6 per bp.

Hudson and Kaplan (1995) show that the effect on a neutral marker at the center of a region of length $R$
of a collection of deleterious loci each having selective disadvantage $sh$ and total mutation rate $U$ 
spread across the region is to reduce diversity by $\exp(-U/(2sh+R))$,
which for a large region is roughly $\exp(-u/r)$, where $u/r$ is the ratio of deleterious mutation rate to recombination rate.
Here we have $U=5$, $R=0.625$, and $sh=.043$, so $U/(2sh+R)=7$.
The first quarter of the chromosome has $U=3.125$ (and $R=0.156$), so $U/(2sh+R)=12.9$,
while the last quarter has $U=0.395$ and so $U/(2sh+R)=1.63$.
But, we want to know if this will affect *spatial structure* --
we want the variation in $N_e$ to affect whether there is strong spatial structure or not.
To do this, we want the neutral model to have little spatial structure, 
but to have stronger spatial structure with lower $N_e$.

Then, we want the universal ancestor to be maybe another 10,000 generations above the start of the simulation.

## Two-dimensional

Here's what we did: a 6x6 grid for 4,000 generations
(to fit within 32G memory):
```
# elapsed: 12:39:49 / kernel: 19.46 / user: 45569.19 / mem: 30613852
OUTBASE="bground_sim_4000gens_6x6"
/usr/bin/time --format='elapsed: %E / kernel: %S / user: %U / mem: %M' ./background-sim.py -T 4000 -N 100 -w 6 -L 25e6 -l 1000 -m 4e-3 -u 5e-3 -r 2.5e-8 -a .23 -b 5.34 -s ${OUTBASE}.selloci \
            -A 10000 -k 1000 -U 1e-7 -o ${OUTBASE}.vcf -t ${OUTBASE}.trees -g ${OUTBASE}.log  &> time_4000gens_6x6_1000samp.log &
```


Timing experiments:
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
# elapsed: 54:06.38 / kernel: 9.33 / user: 3210.98 / mem: 6215856
/usr/bin/time --format='elapsed: %E / kernel: %S / user: %U / mem: %M' ./background-sim.py -T 400 -N 100 -w 8 -L 25e6 -l 1000 -m 4e-3 -u 5e-3 -r 2.5e-8 -a .23 -b 5.34 -s bground_sim_short.selloci \
            -A 10000 -k 1000 -U 1e-7 -o bground_sim_short.vcf -t bground_sim_short.trees -g bground_sim_short.log  &> time_400gens_8x8_1000samp.log
```


With the following:
```
/usr/bin/time --format='elapsed: %E / kernel: %S / user: %U / mem: %M' ./background-sim.py -T 1000 -N 100 -w 8 -L 25e6 -l 1000 -m 4e-3 -u 5e-3 -r 2.5e-8 -a .23 -b 5.34 -s bground_sim_short.selloci \
            -A 10000 -k 1000 -U 1e-7 -o bground_sim_short.vcf -t bground_sim_short.trees -g bground_sim_short.log
```
Memory: 16Gb for 1000 gens.
2h to simulate; 1h to simplify.

Note: "mem" is "Maximum resident set size of the process during its lifetime, in Kilobytes."
```
# 100 gens: 
# elapsed: 12:04.53 / kernel: 14.59 / user: 623.08 / mem: 2864400
# 7 min to sim, 2 min to simplify
NGENS=100
OUTBASE="bground_sim_${NGENS}gens_6x6"
/usr/bin/time --format='elapsed: %E / kernel: %S / user: %U / mem: %M' ./background-sim.py -T $NGENS -N 100 -w 8 -L 25e6 -l 1000 -m 4e-3 -u 5e-3 -r 2.5e-8 -a .23 -b 5.34 -s ${OUTBASE}.selloci \
            -A 10000 -k 1000 -U 1e-7 -o ${OUTBASE}.vcf -t ${OUTBASE}.trees -g ${OUTBASE}.log  &> time_${OUTBASE}.log

# 200 gens: 
# elapsed: 22:41.62 / kernel: 9.42 / user: 1300.16 / mem: 3818808
# 17 min to sim, 3 min to simplify
NGENS=200
OUTBASE="bground_sim_${NGENS}gens_6x6"
/usr/bin/time --format='elapsed: %E / kernel: %S / user: %U / mem: %M' ./background-sim.py -T $NGENS -N 100 -w 8 -L 25e6 -l 1000 -m 4e-3 -u 5e-3 -r 2.5e-8 -a .23 -b 5.34 -s ${OUTBASE}.selloci \
            -A 10000 -k 1000 -U 1e-7 -o ${OUTBASE}.vcf -t ${OUTBASE}.trees -g ${OUTBASE}.log  &> time_${OUTBASE}.log

# 300 gens: 
# elapsed: 35:30.72 / kernel: 9.09 / user: 2103.78 / mem: 5059304
# 29 min to sim, 6 min to simplify
NGENS=300
OUTBASE="bground_sim_${NGENS}gens_6x6"
/usr/bin/time --format='elapsed: %E / kernel: %S / user: %U / mem: %M' ./background-sim.py -T $NGENS -N 100 -w 8 -L 25e6 -l 1000 -m 4e-3 -u 5e-3 -r 2.5e-8 -a .23 -b 5.34 -s ${OUTBASE}.selloci \
            -A 10000 -k 1000 -U 1e-7 -o ${OUTBASE}.vcf -t ${OUTBASE}.trees -g ${OUTBASE}.log  &> time_${OUTBASE}.log

# 400 gens: 
# elapsed: 54:06.38 / kernel: 9.33 / user: 3210.98 / mem: 6215856
# 40 min to sim, 9 min to simplify
NGENS=400
OUTBASE="bground_sim_${NGENS}gens_6x6"
/usr/bin/time --format='elapsed: %E / kernel: %S / user: %U / mem: %M' ./background-sim.py -T $NGENS -N 100 -w 8 -L 25e6 -l 1000 -m 4e-3 -u 5e-3 -r 2.5e-8 -a .23 -b 5.34 -s ${OUTBASE}.selloci \
            -A 10000 -k 1000 -U 1e-7 -o ${OUTBASE}.vcf -t ${OUTBASE}.trees -g ${OUTBASE}.log  &> time_${OUTBASE}.log
```

## One-dimensional

With population size of $N=100$ and migration rate of $m=4e-3$,
coalescence occurs before migration almost surely.
In fact we want little spatial structure in the absence of background selection,
so take $m=0.1$.

To fit within 32G memory, we can do a linear arrangement of 25 populations for for 5,000 generations:

Testing:
```
# elapsed: 4:27.36 / kernel: 23.94 / user: 200.72 / mem: 1162724
OUTBASE="bground_sim_50gens_25x1"
/usr/bin/time --format='elapsed: %E / kernel: %S / user: %U / mem: %M' ./background-sim.py -T 50 -N 100 -w 25 -y 1 -L 25e6 -l 1000 -m 0.1 -u 5e-3 -r 2.5e-8 -a .23 -b 5.34 -s ${OUTBASE}.selloci \
            -A 10000 -k 1000 -U 1e-7 -o ${OUTBASE}.vcf -t ${OUTBASE}.trees -g ${OUTBASE}.log  &> time_${OUTBASE}.log

# elapsed: 6:30.24 / kernel: 11.55 / user: 345.27 / mem: 1174148
OUTBASE="bground_sim_100gens_25x1"
/usr/bin/time --format='elapsed: %E / kernel: %S / user: %U / mem: %M' ./background-sim.py -T 100 -N 100 -w 25 -y 1 -L 25e6 -l 1000 -m 0.1 -u 5e-3 -r 2.5e-8 -a .23 -b 5.34 -s ${OUTBASE}.selloci \
            -A 10000 -k 1000 -U 1e-7 -o ${OUTBASE}.vcf -t ${OUTBASE}.trees -g ${OUTBASE}.log  &> time_${OUTBASE}.log

# elapsed: 14:06:10 / kernel: 29.68 / user: 50726.24 / mem: 26792648
OUTBASE="bground_sim_5000gens_25x1"
OUTDIR=$OUTBASE
mkdir -p $OUTDIR
/usr/bin/time --format='elapsed: %E / kernel: %S / user: %U / mem: %M' ./background-sim.py -T 5000 -N 100 -w 25 -y 1 -L 25e6 -l 1000 -m 0.1 -u 5e-3 -r 2.5e-8 -a .23 -b 5.34 -s ${OUTDIR}/${OUTBASE}.selloci \
            -A 10000 -k 1000 -U 1e-7 -o ${OUTDIR}/${OUTBASE}.vcf -t ${OUTDIR}/${OUTBASE}.trees -g ${OUTDIR}/${OUTBASE}.log  &> ${OUTDIR}/time_${OUTBASE}.log
```

## Other parameters

Higher migration rate, but still two-dimensional:
```
# elapsed: 12:44:59 / kernel: 13.86 / user: 45883.72 / mem: 26796140
OUTBASE="bground_sim_5000gens_5x5"
OUTDIR=$OUTBASE
mkdir -p $OUTDIR
/usr/bin/time --format='elapsed: %E / kernel: %S / user: %U / mem: %M' ./background-sim.py -T 5000 -N 100 -w 5 -y 5 -L 25e6 -l 1000 -m 0.1 -u 5e-3 -r 2.5e-8 -a .23 -b 5.34 -s ${OUTDIR}/${OUTBASE}.selloci \
            -A 10000 -k 1000 -U 1e-7 -o ${OUTDIR}/${OUTBASE}.vcf -t ${OUTDIR}/${OUTBASE}.trees -g ${OUTDIR}/${OUTBASE}.log  &> ${OUTDIR}/time_${OUTBASE}.log
```

Not quite so high migration rate, and higher recombination rate:
```
#
OUTBASE="bground_sim_5000gens_5x5_migr0.01_recomb_1e-7"
OUTDIR=$OUTBASE
mkdir -p $OUTDIR
/usr/bin/time --format='elapsed: %E / kernel: %S / user: %U / mem: %M' ./background-sim.py -T 5000 -N 100 -w 5 -y 5 -L 25e6 -l 1000 -m 0.01 -u 5e-3 -r 1e-7 -a .23 -b 5.34 -s ${OUTDIR}/${OUTBASE}.selloci \
            -A 10000 -k 1000 -U 1e-7 -o ${OUTDIR}/${OUTBASE}.vcf -t ${OUTDIR}/${OUTBASE}.trees -g ${OUTDIR}/${OUTBASE}.log  &> ${OUTDIR}/time_${OUTBASE}.log
```

Even lower migration rate and more recent MRCA:
```
#
OUTBASE="bground_sim_4000gens_3x3_migr0.001_recomb_1e-7_A_1000"
OUTDIR=$OUTBASE
mkdir -p $OUTDIR
/usr/bin/time --format='elapsed: %E / kernel: %S / user: %U / mem: %M' ./background-sim.py -T 4000 -N 100 -w 3 -y 3 -L 25e6 -l 1000 -m 0.001 -u 5e-3 -r 1e-7 -a .23 -b 5.34 -s ${OUTDIR}/${OUTBASE}.selloci \
            -A 1000 -k 1000 -U 1e-7 -o ${OUTDIR}/${OUTBASE}.vcf -t ${OUTDIR}/${OUTBASE}.trees -g ${OUTDIR}/${OUTBASE}.log  &> ${OUTDIR}/time_${OUTBASE}.log
```

## Comparison with/without msprime

```
# elapsed: 6:30.24 / kernel: 11.55 / user: 345.27 / mem: 1174148
OUTBASE="bground_sim_100gens_25x1"
/usr/bin/time --format='elapsed: %E / kernel: %S / user: %U / mem: %M' ./background-sim.py -T 100 -N 100 -w 25 -y 1 -L 25e6 -l 1000 -m 0.1 -u 5e-3 -r 2.5e-8 -a .23 -b 5.34 -s ${OUTBASE}.selloci \
            -A 10000 -k 1000 -U 1e-7 -o ${OUTBASE}.vcf -t ${OUTBASE}.trees -g ${OUTBASE}.log  &> time_${OUTBASE}.log

# elapsed: 2:27.53 / kernel: 0.56 / user: 147.28 / mem: 650148
OUTBASE="bground_sim_100gens_25x1_nomsprime"
/usr/bin/time --format='elapsed: %E / kernel: %S / user: %U / mem: %M' ./background-sim-no-msprime.py -T 100 -N 100 -w 25 -y 1 -L 25e6 -l 1000 -m 0.1 -u 5e-3 -r 2.5e-8 -a .23 -b 5.34 -s ${OUTBASE}.selloci \
            -A 10000 -k 1000 -U 1e-7 -o ${OUTBASE}.vcf -t ${OUTBASE}.trees -g ${OUTBASE}.log  &> time_${OUTBASE}.log
```

## Population history

### Split-and-merge

Suppose that a set of four populations of diploid size N 
transitions from an (AB)(CD) split to a (AC)(BD) split T generations ago.
Let $m>1/4N$ be the migration rate between 'nearby' populations 
and $M<1/4N$ be the rate between 'separated' populations.
If $\exp(-T/4N)$, the probability of no coalescence since the split, 
is small, then most trees follow (AC)(BD), and conversely.
We want to take $T/4N$ to be smallish, 
but big enough that background selection can make $T/4N_e$ big.

To do this in msprime, see:
```
CHRLEN=5e6
OUTDIR=sim_${CHRLEN}_${RANDOM}
time python3 neutral-split-sim.py -o $OUTDIR --Ne 1000 20000 --nsamples 100 --chrom_length $CHRLEN -T 0.1 --relative_fast_m 1
for vcf in ${OUTDIR}/*.vcf; do 
    CHROM=${vcf%%.vcf}
    CHROM=${CHROM: -1}
    sed -e "s/ID=1,/ID=${CHROM},/" -i $vcf
    sed -e "s/^1/${CHROM}/" -i $vcf
    bcftools convert -O b $vcf -o ${vcf%%vcf}bcf
    bcftools index ${vcf%%vcf}bcf
    rm $vcf
done
printf -v CHRLENBP "%.f" "$CHRLEN"
WINLEN=$((CHRLENBP/100))
for NPC in 1 2 3; do
    (LODIR=${OUTDIR}/bp_${WINLEN}_npc_${NPC};
    ./run_lostruct.R -i ${OUTDIR} -k $NPC -t bp -s $WINLEN -o $LODIR -I ${OUTDIR}/samples.tsv;
    Rscript -e "templater::render_template('summarize_run.Rmd',output='${LODIR}/run-summary.html',change.rootdir=TRUE)")&
done
```

### Three populations

Alternatively, suppose that in recent times three populations with strong asymmetric migration `A<-B|C`
that previously was `A|B->C`.

To do this in msprime (i.e., neutrally, but with different Ne on different chromsomes):
```
CHRLEN=1e6
OUTDIR=threesim_${CHRLEN}_${RANDOM}
time python3 threeway-neutral-split-sim.py -o $OUTDIR --Ne 10000 100000 --nsamples 100 --chrom_length $CHRLEN -T 0.25 --relative_fast_m 1
for vcf in ${OUTDIR}/*.vcf; do 
    CHROM=${vcf%%.vcf}
    CHROM=${CHROM: -1}
    sed -e "s/ID=1,/ID=${CHROM},/" -i $vcf
    sed -e "s/^1/${CHROM}/" -i $vcf
    bcftools convert -O b $vcf -o ${vcf%%vcf}bcf
    bcftools index ${vcf%%vcf}bcf
    rm $vcf
done
printf -v CHRLENBP "%.f" "$CHRLEN"
WINLENBP=$((CHRLENBP/100))
SNPNUM0=$(bcftools stats $OUTDIR/chrom0.bcf | grep "number of SNPs" | head -n 1 | cut -f 4)
SNPNUM1=$(bcftools stats $OUTDIR/chrom1.bcf | grep "number of SNPs" | head -n 1 | cut -f 4)
SNPNUM=$((SNPNUM0+SNPNUM1))
WINLENSNP=$((SNPNUM/200))
for NPC in 1 2 3; do
    # in BP
    (LODIR=${OUTDIR}/bp_${WINLENBP}_npc_${NPC};
    ./run_lostruct.R -i ${OUTDIR} -k $NPC -t bp -s $WINLENBP -o $LODIR -I ${OUTDIR}/samples.tsv;
    Rscript -e "templater::render_template('summarize_run.Rmd',output='${LODIR}/run-summary.html',change.rootdir=TRUE)")&
    # in SNPs
    (LODIR=${OUTDIR}/snp_${WINLENSNP}_npc_${NPC};
    ./run_lostruct.R -i ${OUTDIR} -k $NPC -t snp -s $WINLENSNP -o $LODIR -I ${OUTDIR}/samples.tsv;
    Rscript -e "templater::render_template('summarize_run.Rmd',output='${LODIR}/run-summary.html',change.rootdir=TRUE)")&
done
```

/usr/bin/time --format='elapsed: %E / kernel: %S / user: %U / mem: %M' ./background-sim.py -T 5000 -N 100 -w 25 -y 1 -L 25e6 -l 1000 -m 0.1 -u 5e-3 -r 2.5e-8 -a .23 -b 5.34 -s ${OUTDIR}/${OUTBASE}.selloci \
            -A 10000 -k 1000 -U 1e-7 -o ${OUTDIR}/${OUTBASE}.vcf -t ${OUTDIR}/${OUTBASE}.trees -g ${OUTDIR}/${OUTBASE}.log  &> ${OUTDIR}/time_${OUTBASE}.log

### Gamma-distributed s

To do this in simuPOP (i.e., for real):
```

CHRLEN=25e6
POPSIZE=1000
NLOCI=4000
NSAMPLES=200
OUTDIR=threesim_background_${CHRLEN}_${POPSIZE}_${NLOCI}_${NSAMPLES}_${RANDOM}
OUTBASE=threesim
mkdir -p $OUTDIR
/usr/bin/time --format='elapsed: %E / kernel: %S / user: %U / mem: %M' python3 threeway-background-sim.py \
    -o ${OUTDIR}/${OUTBASE}.vcf -t ${OUTDIR}/${OUTBASE}.trees -g ${OUTDIR}/${OUTBASE}.log \
    --nloci $NLOCI --popsize $POPSIZE --nsamples $NSAMPLES --length $CHRLEN --relative_switch_time 0.1 -T 100 -A 100 \
    --recomb_rate 1e-7 --sel_mut_rate 1e-3 --relative_fast_M 1  --relative_slow_m .01  &> ${OUTDIR}/time_${OUTBASE}.log

/usr/bin/time --format='elapsed: %E / kernel: %S / user: %U / mem: %M' python3 tree-stats.py -t ${OUTDIR}/${OUTBASE}.trees \
    -s ${OUTDIR}/samples.tsv -n 100 -o ${OUTDIR}/divergences.tsv &> ${OUTDIR}/time_divergences.log

printf -v CHRLENBP "%.f" "$CHRLEN"
bcftools convert -O b -o ${OUTDIR}/${OUTBASE}.bcf ${OUTDIR}/${OUTBASE}.vcf
bcftools index ${OUTDIR}/${OUTBASE}.bcf
WINLENBP=$((CHRLENBP/100))
SNPNUM=$(bcftools stats $OUTDIR/${OUTBASE}.bcf | grep "number of SNPs" | head -n 1 | cut -f 4)
WINLENSNP=$((SNPNUM/200))
for NPC in 1 2 3; do
    # in BP
    ( LODIR=${OUTDIR}/bp_${WINLENBP}_npc_${NPC};
      ./run_lostruct.R -i ${OUTDIR} -k $NPC -t bp -s $WINLENBP -o $LODIR -I ${OUTDIR}/samples.tsv;
      Rscript -e "templater::render_template('summarize_run.Rmd',output='${LODIR}/run-summary.html',change.rootdir=TRUE)")&
    # in SNPs
    ( LODIR=${OUTDIR}/snp_${WINLENSNP}_npc_${NPC};
      ./run_lostruct.R -i ${OUTDIR} -k $NPC -t snp -s $WINLENSNP -o $LODIR -I ${OUTDIR}/samples.tsv;
      Rscript -e "templater::render_template('summarize_run.Rmd',output='${LODIR}/run-summary.html',change.rootdir=TRUE)")&
done
    
```

This one orignally ran for `(0.25*3*2*1000 = 1500) + 100` generations with 3000 individuals
modified it to run for `(0.25*3*2*1200 + 1800) = 3600` generations with `3*1200=3600` individuals
```
RUNID=$RANDOM
OUTDIR=threesim_background_$RUNID
mkdir -p $OUTDIR
/usr/bin/time --format='elapsed: %E / kernel: %S / user: %U / mem: %M' python3 threeway-background-sim.py \
   --ancestor_age 100.0  \
    --gamma_alpha 0.23  \
    --gamma_beta 5.34  \
    --generations 1800  \
    --length 25000000  \
    --mut_rate 1e-07  \
    --nloci 1000  \
    --nsamples 500  \
    --popsize 1200  \
    --recomb_rate 2.5e-8  \
   --relative_fast_M 1.0  \
    --relative_slow_m 0.1  \
    --relative_switch_time 0.25  \
    --sel_mut_rate 5e-3  \
   --logfile "$OUTDIR/threesim.log"  \
   --outfile "$OUTDIR/threesim.vcf"  \
   --selloci_file "$OUTDIR/threesim.selloci"  \
   --treefile "$OUTDIR/threesim.trees" &> $OUTDIR/time_threesim.log
```
and this one is similar but smaller slow_m:
```
/usr/bin/time --format='elapsed: %E / kernel: %S / user: %U / mem: %M' python3 threeway-background-sim.py \
    --ancestor_age 100.0 
    --gamma_alpha 0.23 
    --gamma_beta 5.34 
    --generations 100 \
    --length 25000000.0 
    --mut_rate 1e-07 
    --nloci 4000 
    --nsamples 200 \
    --popsize 1000 \
    --recomb_rate 1e-07 
    --relative_fast_M 1.0 
    --relative_slow_m 0.01 
    --relative_switch_time 0.25 \
    --logfile 'threesim_background_25e6_1000_4000_200_22290/threesim.log' \
    --outfile 'threesim_background_25e6_1000_4000_200_22290/threesim.vcf' 
    --sel_mut_rate 0.001 
    --samples_file 'threesim_background_25e6_1000_4000_200_22290/samples.tsv' 
    --selloci_file 'threesim_background_25e6_1000_4000_200_22290/sel_loci.txt' \
    --treefile 'threesim_background_25e6_1000_4000_200_22290/threesim.trees'

    --fast_M 0.00025 
    --slow_m 2.5e-06 
    --switch_time 1500 \

```

### Fixed *s*

```

CHRLEN=0.5e7
POPSIZE=500
NLOCI=4000
SELCOEF=0.1
OUTDIR=threesim_bg_fixed_s_${CHRLEN}_${POPSIZE}_${NLOCI}_${SELCOEF}_${RANDOM}
NSAMPLES=100
OUTBASE=threesim
mkdir -p $OUTDIR
/usr/bin/time --format='elapsed: %E / kernel: %S / user: %U / mem: %M' python3 threeway-background-fixed-s-sim.py \
    -o ${OUTDIR}/${OUTBASE}.vcf -t ${OUTDIR}/${OUTBASE}.trees -g ${OUTDIR}/${OUTBASE}.log \
    --nloci $NLOCI --popsize $POPSIZE --nsamples $NSAMPLES --length $CHRLEN --relative_switch_time 0.1 -T 100 -A 100 --mut_rate 1e-6 \
    --selection_coef $SELCOEF --recomb_rate 1e-7 --sel_mut_rate 1e-3 --relative_fast_M 1  --relative_slow_m .01  &> ${OUTDIR}/time_${OUTBASE}.log

/usr/bin/time --format='elapsed: %E / kernel: %S / user: %U / mem: %M' python3 tree-stats.py -t ${OUTDIR}/${OUTBASE}.trees \
    -s ${OUTDIR}/samples.tsv -n 100 -o ${OUTDIR}/divergences.tsv &> ${OUTDIR}/time_divergences.log

printf -v CHRLENBP "%.f" "$CHRLEN"
bcftools convert -O b -o ${OUTDIR}/${OUTBASE}.bcf ${OUTDIR}/${OUTBASE}.vcf
bcftools index ${OUTDIR}/${OUTBASE}.bcf
WINLENBP=$((CHRLENBP/100))
SNPNUM=$(bcftools stats $OUTDIR/${OUTBASE}.bcf | grep "number of SNPs" | head -n 1 | cut -f 4)
WINLENSNP=$((SNPNUM/200))
for NPC in 1 2 3; do
    # in BP
    ( LODIR=${OUTDIR}/bp_${WINLENBP}_npc_${NPC};
      ./run_lostruct.R -i ${OUTDIR} -k $NPC -t bp -s $WINLENBP -o $LODIR -I ${OUTDIR}/samples.tsv;
      Rscript -e "templater::render_template('summarize_run.Rmd',output='${LODIR}/run-summary.html',change.rootdir=TRUE)")&
    # in SNPs
    ( LODIR=${OUTDIR}/snp_${WINLENSNP}_npc_${NPC};
      ./run_lostruct.R -i ${OUTDIR} -k $NPC -t snp -s $WINLENSNP -o $LODIR -I ${OUTDIR}/samples.tsv;
      Rscript -e "templater::render_template('summarize_run.Rmd',output='${LODIR}/run-summary.html',change.rootdir=TRUE)")&
done
    
```



# Testing simuPOP

We'll try simulating a basic popultaion with evenly distributed loci and a range of fixed selection coefficients.

As in Hudson & Kaplan 1995, with

- $U$ = total deleterious mutation rate (= number of loci times per-locus rate)
- $sh$ = effect of a single mutation,

we have that:

- the number of deleterious mutations per individual is Poisson($U/(2sh)$)
- without recombination, divergence is $\pi_0 \exp(-U/(2sh))$ (Charlesworth et al 1993)
- with recombination distance between the focal locus and locus $x$, 
    divergence is $\pi_0 \exp\left(-\int u(x) sh dx / 2(sh + R(x)\right)^2$;
    which in the center of the chromsome with constant recombination and mutation
    is $\pi_0 \exp\left(-U/(2sh+R)\right)$, where $R$ is the length of the chromosome in Morgans

Below, $U$=`NSEL * SELMUTRATE` and $sh$=`SELCOEF` and $R=10^6$` * RECOMBRATE` and $\pi_0$=`2 POPSIZE`.

Trial baloon:
```sh

runsim () {
    POPSIZE=$1
    NSEL=$2
    SELCOEF=$3
    SELMUTRATE=$4
    RECOMBRATE=$5
    OUTDIR="test_${POPSIZE}_${NSEL}_${SELCOEF}_${SELMUTRATE}_${RECOMBRATE}_${RANDOM}"; mkdir -p $OUTDIR
    /usr/bin/time --format='elapsed: %E / kernel: %S / user: %U / mem: %M' \
        python3 simple-background-sim.py --generations $((5 * $POPSIZE)) --popsize $POPSIZE --length 1e6 --nloci $NSEL --recomb_rate $RECOMBRATE \
            --sel_mut_rate $SELMUTRATE --selection_coef $SELCOEF \
            --nsamples 20 --ancestor_age 100 --mut_rate 0 --treefile $OUTDIR/sim.trees --logfile $OUTDIR/sim.log \
            &> $OUTDIR/time.log
    python3 ../tree-stats.py --treefile $OUTDIR/sim.trees --samples_file $OUTDIR/samples.tsv --n_window 100 --outfile $OUTDIR/divergences.tsv
    echo $OUTDIR
}

runsim 100 100 0.1 1e-3 1e-5
runsim 100 100 0.01 1e-3 1e-5
runsim 100 100 0.1 1e-4 1e-5

runsim 200 400 0.01 1e-3 1e-5
runsim 200 400 0.0 1e-3 1e-5
runsim 200 400 0.01 1e-4 1e-5

runsim 100 1000 0.01 1e-3 1e-5
runsim 100 1000 0.0 1e-3 1e-5
runsim 100 1000 0.01 1e-4 1e-5
runsim 100 1000 0.001 1e-3 1e-5
runsim 100 1000 0.001 1e-4 1e-5

runsim 100 10000 0.01 1e-3 1e-5
runsim 100 10000 0.0 1e-3 1e-5
runsim 100 10000 0.01 1e-4 1e-5
runsim 100 10000 0.001 1e-3 1e-5
runsim 100 10000 0.001 1e-4 1e-5

runsim 200 4000 0.1 1e-3 1e-5
runsim 200 4000 0.01 1e-3 1e-5
runsim 200 4000 0.001 1e-3 1e-5
runsim 200 4000 0.01 1e-4 1e-5

# example had 2382 muts in first 5MB of sequence (=0.5 M)
# runsim 1000 2400 0.01 .001 0.5e-6 &
runsim 500 2400 0.01 .001 0.5e-6 
runsim 500 2400 0.1 .001 0.5e-6 
# and 281 in the last 5MB
# runsim 1000 280 0.01 .001 0.5e-6 &
runsim 500 280 0.01 .001 0.5e-6 &
runsim 500 280 0.1 .001 0.5e-6 &
# and 4000 across 25MB (=2.5M)
#   ... which is an average of 800 per 0.5M
# runsim 1000 800 0.01 .001 0.5e-6 &
runsim 500 800 0.01 .001 0.5e-6 &
runsim 500 800 0.1 .001 0.5e-6 &

# ... what if we did similar on one-tenth the sequence?
runsim 500 800 0.01 .001 0.05e-6 
runsim 500 800 0.1 .001 0.05e-6 

# print name, average divergence
( echo "name divergence popsize U denom length s pi N/div div/theory"; 
    ( for x in test_*; do echo $x $(if [ -f $x/divergences.tsv ]; then cat $x/divergences.tsv | tail -n +2  | awk 'BEGIN { x=0; n=0 } {x+=$3; n+=1} END { print x/n}' 2>/dev/null; else echo "xxx"; fi); done ) | awk -F "_" '{ N=$2; U=$3*$5; R=$6; D=(2*$4+1000000*R); print $0,N,U,D,1000000*R,$4,4*$2*exp(-U/D)}' | awk '{print $0,4*$3/$2,$2/$8}' | sort  -k 9 -n ) | column -t > results.tsv

```


With uneven loci:
```

runsim2 () {
    POPSIZE=$1
    NSEL=$2
    SELCOEF=$3
    SELMUTRATE=$4
    RECOMBRATE=$5
    OUTDIR="uneven_${NSEL}_${SELCOEF}_${SELMUTRATE}_${RECOMBRATE}_${RANDOM}"; mkdir -p $OUTDIR
    /usr/bin/time --format='elapsed: %E / kernel: %S / user: %U / mem: %M' \
        python3 simple-uneven-loci-background-sim.py --generations $((5 * $POPSIZE)) --popsize $POPSIZE --length 1e6 --nloci $NSEL --recomb_rate $RECOMBRATE \
            --sel_mut_rate $SELMUTRATE --selection_coef $SELCOEF \
            --nsamples 20 --ancestor_age 100 --mut_rate 0 --treefile $OUTDIR/sim.trees --logfile $OUTDIR/sim.log \
            &> $OUTDIR/time.log
    python3 ../tree-stats.py --treefile $OUTDIR/sim.trees --samples_file $OUTDIR/samples.tsv --n_window 100 --outfile $OUTDIR/divergences.tsv
    echo $OUTDIR
}

runsim2 100 100 0.1 1e-3 1e-5 & 

# example had 2382 muts in first 5MB of sequence (=0.5 M)
runsim2 500 2400 0.01 .001 0.5e-6 &
runsim2 500 2400 0.1 .001 0.5e-6 &
# and 281 in the last 5MB
runsim2 500 280 0.01 .001 0.5e-6 &
runsim2 500 280 0.1 .001 0.5e-6 &
# and 4000 across 25MB (=2.5M)
#   ... which is an average of 800 per 0.5M
runsim2 500 800 0.01 .001 0.5e-6 &
runsim2 500 800 0.1 .001 0.5e-6 &

# print name, average divergence
( echo "name divergence popsize U denom length s pi N/div div/theory"; 
    ( for x in uneven_*; do echo $x $(if [ -f $x/divergences.tsv ]; then cat $x/divergences.tsv | tail -n +2  | awk 'BEGIN { x=0; n=0 } {x+=$3; n+=1} END { print (n+0==0)?"NA":(x/n)}' 2>/dev/null; else echo "xxx"; fi); done ) | awk -F "_" '{ N=$2; U=$3*$5; R=$6; D=(2*$4+1000000*R); print $0,N,U,D,1000000*R,$4,(D+0==0)?"NA":(4*$2*exp(-U/D))}' | awk '{X=($2+0==0)?"NA":(4*$3/$2); Y=($8+0==0)?"NA":($2/$8); print $0,X,Y}' | sort  -k 9 -n ) | column -t > results_uneven.tsv

```

# Local adaptation

```

localsim () {
    POPSIZE=$1
    NSEL=$2
    SELECTION_COEF=$3
    RELATIVE_M=$4
    RECOMB_RATE=$5
    SEED=$RANDOM
    SCRIPT=local-fixed-s-sim.py 
    OUTDIR="local_${1}_${2}_${3}_${4}_${5}_${SEED}"; mkdir -p $OUTDIR
    /usr/bin/time --format='elapsed: %E / kernel: %S / user: %U / mem: %M' \
        python3 $SCRIPT \
                     --relative_m $RELATIVE_M \
                     --generations $((5 * $POPSIZE)) \
                     --popsize $POPSIZE \
                     --length 1e6  \
                     --nloci $NSEL \
                     --sel_mut_rate 1e-3 \
                     --recomb_rate $RECOMB_RATE \
                     --selection_coef $SELECTION_COEF \
                     --nsamples 20 \
                     --ancestor_age 100 \
                     --mut_rate 1e-5  \
                     --seed $SEED \
                     --treefile $OUTDIR/sim.trees  \
                     --outfile $OUTDIR/sim.vcf \
                     --logfile $OUTDIR/sim.log  \
            &> $OUTDIR/time.log
    python3 ../tree-stats.py --treefile $OUTDIR/sim.trees --samples_file $OUTDIR/samples.tsv --n_window 100 --outfile $OUTDIR/divergences.tsv
    echo $OUTDIR
}

# testing:
localsim 20 10 0.001 1.0 .000001

# elapsed: 0:41.14 / kernel: 0.65 / user: 40.70 / mem: 410508
localsim 100 1000 0.001 1.0 .000001

# elapsed: 2:21.83 / kernel: 0.86 / user: 141.49 / mem: 1113744
localsim 100 1000 0.001 1.0 .00001

# elapsed: 4:34.72 / kernel: 0.92 / user: 274.30 / mem: 1493748
localsim 200 1000 0.001 1.0 .000001

# elapsed: 9:25.62 / kernel: 1.29 / user: 564.74 / mem: 2423340
localsim 200 1000 0.001 1.0 .000004


```


Here there are three populations, connected by migration rates $m = m_{12} < m_{23} = M$;
with selection coefficients opposite in the third population.

- On the neutral half, the migration rates will be as given; if $M > 1/N$ 
  then divergence times between 2 and 3 will be like $1/N$, and if $m < 1/N$
  then divergence times between 1 and 2 or 3 will be like $1/N + 1/m$.

- On the selected half, with $L$ loci of selection coefficient $s$,
  of which half are positive and half are negative,
  a migrant has disadvantage around $Ls/2$, thus replacing $m_{23}$ by $M \exp(-LS/2)$.
  If this is much less than $m$, then on the selected half, 1 and 2 should be closer than to 3.

- 

```

localthree () {
    POPSIZE=$1
    NSEL=$2
    SELECTION_COEF=$3
    SLOW_M=$4
    FAST_M=$5
    RECOMB_RATE=$6
    RELATIVE_GENS=$7
    SEED=$RANDOM
    SCRIPT=threeway-local-fixed-s-sim.py
    OUTDIR="threeway_chroms_${1}_${2}_${3}_${4}_${5}_${6}_${7}_${SEED}"; mkdir -p $OUTDIR
    /usr/bin/time --format='elapsed: %E / kernel: %S / user: %U / mem: %M' \
        python3 $SCRIPT \
                     --relative_m $SLOW_M \
                     --relative_M $FAST_M \
                     --generations $(($RELATIVE_GENS * $POPSIZE)) \
                     --popsize $POPSIZE \
                     --length 1e6  \
                     --nloci $NSEL \
                     --sel_mut_rate 1e-3 \
                     --recomb_rate $RECOMB_RATE \
                     --selection_coef $SELECTION_COEF \
                     --nsamples 20 \
                     --ancestor_age 100 \
                     --mut_rate 1e-5  \
                     --seed $SEED \
                     --treefile $OUTDIR/sim.trees  \
                     --outfile $OUTDIR/sim.vcf \
                     --logfile $OUTDIR/sim.log  \
            &> $OUTDIR/time.log
    echo "Now computing tree stats." >> $OUTDIR/time.log
    /usr/bin/time --format='elapsed: %E / kernel: %S / user: %U / mem: %M' \
        python3 ../tree-stats.py --treefile $OUTDIR/sim.trees --samples_file $OUTDIR/samples.tsv \
        --n_window 100 --outfile $OUTDIR/divergences.tsv &>> $OUTDIR/time.log
    echo $OUTDIR
}

# elapsed: 15:02.56 / kernel: 1.72 / user: 901.34 / mem: 3587288
localthree 200 1000 0.001 1.0 0.1 .000004 10 &

localthree 200 1000 0.001 5.0 0.02 .000004 10 &

localthree 200 1000 0.003 5.0 0.2 .000004 10 &

localthree 200 1000 0.003 5.0 0.5 .000004 10 &

# increase recomb, lower selection
localthree 100 100 0.003 5.0 0.5 .000004 10 &

# fewer selected sites?
localthree 100 100 0.003 5.0 0.5 .000004 10 &

# increase migration rate?
localthree 100 500 0.006 20.0 2.0 .000004 10 &
localthree 100 500 0.006 20.0 5.0 .000004 10 &

```

## Here is a more symmetric situation:

```

symthree () {
    POPSIZE=$1
    NSEL=$2
    SELECTION_COEF=$3
    RELATIVE_M=$4
    RECOMB_RATE=$5
    RELATIVE_GENS=$6
    SEED=$RANDOM
    SCRIPT=threeway-symmetric-local-fixed-s-sim.py
    OUTDIR="threeway_sym_${1}_${2}_${3}_${4}_${5}_${6}_${SEED}"; mkdir -p $OUTDIR
    /usr/bin/time --format='elapsed: %E / kernel: %S / user: %U / mem: %M' \
        python3 $SCRIPT \
                     --relative_m $RELATIVE_M \
                     --generations $(($RELATIVE_GENS * $POPSIZE)) \
                     --popsize $POPSIZE \
                     --length 1e6  \
                     --nloci $NSEL \
                     --sel_mut_rate 1e-3 \
                     --recomb_rate $RECOMB_RATE \
                     --selection_coef $SELECTION_COEF \
                     --nsamples 20 \
                     --ancestor_age 100 \
                     --mut_rate 1e-5  \
                     --seed $SEED \
                     --treefile $OUTDIR/sim.trees  \
                     --outfile $OUTDIR/sim.vcf \
                     --logfile $OUTDIR/sim.log  \
            &> $OUTDIR/time.log
    echo "Now computing tree stats." >> $OUTDIR/time.log
    /usr/bin/time --format='elapsed: %E / kernel: %S / user: %U / mem: %M' \
        python3 ../tree-stats.py --treefile $OUTDIR/sim.trees --samples_file $OUTDIR/samples.tsv \
        --n_window 100 --outfile $OUTDIR/divergences.tsv &>> $OUTDIR/time.log
    echo $OUTDIR
}

# testing
symthree 20 100 0.001 1.0 .000004 10 &

# 
symthree 200 1000 0.001 10.0 .000004 10 &

# larger Ns
symthree 200 1000 0.005 10.0 .000004 10 &
symthree 200 1000 0.010 10.0 .000004 10 &

```


## Looking at results:
```
divs <- file.path(list.files(".", "threeway", full.names=TRUE), "divergences.tsv")
divs <- divs[file.exists(divs)]
divs <- divs[order(file.info(divs)$ctime)][-(1:2)]
x <- lapply(divs, read.table, header=TRUE)
names(x) <- basename(dirname(divs))
cols <- c(grey(.8), 'red', 'green', 'red', grey(.5), 'purple', 'green', 'purple', grey(.3))
names(cols) <- paste0("X", c("0_0", "1_0", "2_0", "0_1", "1_1", "2_1", "0_2", "1_2", "2_2"))
layout(matrix(1:24,nrow=4))
for (k in seq_along(x)) { 
    xx <- x[[k]]
    matplot((xx[,1]+xx[,2])/2,xx[,-(1:2)], type='l', lty=1, col=cols[colnames(xx)[-(1:2)]],
        main=names(x)[k], xlab='position (bp)') 
    legend("topright", lty=1, col=cols[colnames(xx)[-(1:2)]], legend=colnames(xx)[-(1:2)])
}

```

lostructify:

```

lostruct () {
    OUTDIR=$1
    WINLENBP=$((1000000/40))
    bcftools convert -O b -o ${OUTDIR}/sim.bcf ${OUTDIR}/sim.vcf
    bcftools index ${OUTDIR}/sim.bcf
    for NPC in 1 2 3; do
        (LODIR=${OUTDIR}/bp_${WINLENBP}_npc_${NPC};
        ./run_lostruct.R -i ${OUTDIR} -k $NPC -t bp -s $WINLENBP -o $LODIR -I ${OUTDIR}/samples.tsv;
        Rscript -e "templater::render_template('summarize_run.Rmd',output='${LODIR}/run-summary.html',change.rootdir=TRUE)")&
    done
}

```
