#!/usr/bin/env python3
description = '''
Simulates and writes to msprime/vcf, three populations connected by migration,
with a number of SNPs under selection evenly distributed on the chromosome.
The left two pops share an environment that the left half the chromosome affects, 
so the right pop has opposite selection coefficients,
while the right two share an environment affected by SNPs on the right half.
'''

import gzip
import sys, os
import math
import time
import random
from ftprime import RecombCollector, ind_to_chrom, mapa_labels
import msprime
import argparse

def fileopt(fname,opts):
    '''Return the file referred to by fname, open with options opts;
    if fname is "-" return stdin/stdout; if fname ends with .gz run it through gzip.
    '''
    if fname == "-":
        if opts == "r":
            fobj = sys.stdin
        elif opts == "w":
            fobj = sys.stdout
        else:
            print("Something not right here.")
    elif fname[len(fname)-3:len(fname)]==".gz":
        fobj = gzip.open(fname,opts)
    else:
        fobj = open(fname,opts)
    return fobj

parser = argparse.ArgumentParser(description=description)
parser.add_argument('--relative_m', '-m', default=1.0, type=float, 
        help="Migration rate between pops in units of population size.")
parser.add_argument("--generations", "-T", type=int, help="number of generations to run for.")
parser.add_argument("--popsize", "-N", type=int, help="size of each subpopulation", default=100)
parser.add_argument("--length", "-L", type=float, help="number of bp in the chromosome", default=1e4)
parser.add_argument("--nloci", "-l", type=int, help="number of selected loci", default=20)
parser.add_argument("--sel_mut_rate", "-u", type=float, help="mutation rate of selected alleles", default=1e-7)
parser.add_argument("--recomb_rate", "-r", type=float, help="recombination rate", default=1e-7)
parser.add_argument("--selection_coef","-S", type=float, help="strength of selection",default=.1)
parser.add_argument("--nsamples", "-k", type=int, help="number of *diploid* samples, total")
parser.add_argument("--ancestor_age", "-A", type=float, help="time to ancestor above beginning of sim")
parser.add_argument("--mut_rate", "-U", type=float, help="mutation rate", default=1e-7)
parser.add_argument("--seed", "-d", type=int, help="random seed", default=random.randrange(1,1000))

parser.add_argument("--treefile", "-t", help="name of output file for trees (default: not output)", default=None)
parser.add_argument("--outfile", "-o", help="name of output VCF file (default: not output)", default=None)
parser.add_argument("--logfile", "-g", help="name of log file (or '-' for stdout)", default="-")
parser.add_argument("--selloci_file", "-s", help="name of file to output selected locus information (default: (dir)/sel_loci.txt)")
parser.add_argument("--samples_file", "-e", help="name of file to output information on samples (default=(dir)/samples.tsv)")

args = parser.parse_args()

import simuOpt
import simuPOP as sim

sim.setRNG(seed=args.seed)
random.seed(args.seed)

if args.outfile is not None:
    outfile = fileopt(args.outfile, "w")
logfile = fileopt(args.logfile, "w")
if args.selloci_file is None:
    args.selloci_file = os.path.join(os.path.dirname(args.treefile),"sel_loci.txt")
selloci_file = args.selloci_file
if args.samples_file is None:
    args.samples_file = os.path.join(os.path.dirname(args.treefile),"samples.tsv")
samples_file = fileopt(args.samples_file,"w")

# compute these here so they get recorded in the log
args.m = args.relative_m/(4*args.popsize)

logfile.write("Options:\n")
logfile.write(str(args)+"\n")
logfile.write(time.strftime('%X %x %Z')+"\n")
logfile.write("----------\n")
logfile.flush()

npops=3

# loci evenly spaced on chromosome
# except must have one on each end
locus_position=[x*args.length/args.nloci for x in range(args.nloci-1)] + [args.length]

# initially polymorphic alleles
init_freqs=[[k/100,1-k/100,0,0] for k in range(10)]
locus_classes=[min(len(init_freqs)-1,math.floor(random.expovariate(1))) for k in range(args.nloci)]
init_classes=[list(filter(lambda k: locus_classes[k]==x,range(args.nloci))) for x in range(len(init_freqs))]

logfile.write("Locus positions:\n")
logfile.write(str(locus_position)+"\n")
logfile.write("----------\n")
logfile.flush()

init_geno=[sim.InitGenotype(freq=init_freqs[k],loci=init_classes[k]) for k in range(len(init_freqs))]

class PlusMinusFitness:
    '''
    left half the chromosome is governed by two regimes: A/a
    right half the chromosome is governed by: B/b
    spatial arrangement is:
    A A a
    B b b
    1 2 3
    We do this by flipping signs for lower-case letters.
    '''
    def __init__(self, s, a_cutoff, b_cutoff):
        self.coefMap = {}
        self.s = s
        self.a_cutoff = a_cutoff
        self.b_cutoff = b_cutoff
        self.AB_subpops = [0]
        self.Ab_subpops = [1]
        self.ab_subpops = [2]
    def AB(self, loc, alleles):
        # because sign is assigned for each locus, we need to make sure the
        # same sign is used for fitness of genotypes 01 (1-s) and 11 (1-2s)
        # at each locus
        if loc in self.coefMap:
            s = self.coefMap[loc]
        else:
            s = random.choice([-1.0,1.0]) * self.s
            self.coefMap[loc] = s
        # print(str(loc)+":"+str(alleles)+"\n")
        # needn't return fitness for alleles=(0,0) as simupop knows that's 1
        if 0 in alleles:
            return max(0.0, 1. - s)
        else:
            return max(0.0, 1. - 2.*s)
    def Ab(self, loc, alleles):
        if loc in self.coefMap:
            s = self.coefMap[loc]
        else:
            s = random.choice([-1.0,1.0]) * self.s
            self.coefMap[loc] = s
        if loc > self.b_cutoff:
            s = (-1) * s
        if 0 in alleles:
            return max(0.0, 1. - s)
        else:
            return max(0.0, 1. - 2.*s)
    def ab(self, loc, alleles):
        if loc in self.coefMap:
            s = self.coefMap[loc]
        else:
            s = random.choice([-1.0,1.0]) * self.s
            self.coefMap[loc] = s
        if loc > self.b_cutoff:
            s = (-1) * s
        if loc > self.a_cutoff:
            s = (-1) * s
        if 0 in alleles:
            return max(0.0, 1. - s)
        else:
            return max(0.0, 1. - 2.*s)


fitness = PlusMinusFitness(args.selection_coef, args.nloci/2, arg.nloci/2)

pop = sim.Population(
        size=[args.popsize]*npops, 
        loci=[args.nloci,2], 
        lociPos=locus_position,
        infoFields=['ind_id','fitness','migrate_to'])

id_tagger = sim.IdTagger()
id_tagger.apply(pop)

# record recombinations
rc = RecombCollector(
        first_gen=pop.indInfo("ind_id"), ancestor_age=args.ancestor_age, 
                              length=2*args.length, 
                              locus_position=locus_position)

migr_mat = [[ 0, args.m, 0 ],
            [ args.m, 0, args.m ],
            [ 0, args.m, 0 ] ]

pop.evolve(
    initOps=[
        sim.InitSex(),
    ]+init_geno,
    preOps=[
        sim.PyOperator(lambda pop: rc.increment_time() or True),
        sim.Migrator(
            rate=migr_mat,
            mode=sim.BY_PROBABILITY),
        sim.SNPMutator(u=args.sel_mut_rate, v=args.sel_mut_rate),
        sim.PyMlSelector(fitness.AB, subPops=fitness.AB_subpops,
            output=">>"+selloci_file),
        sim.PyMlSelector(fitness.Ab, subPops=fitness.Ab_subpops,
            output=">>"+selloci_file),
        sim.PyMlSelector(fitness.ab, subPops=fitness.ab_subpops,
            output=">>"+selloci_file),
    ],
    matingScheme=sim.RandomMating(
        ops=[
            id_tagger,
            sim.Recombinator(intensity=args.recomb_rate,
                output=rc.collect_recombs,
                infoFields="ind_id"),
        ] ),
    postOps=[
        sim.Stat(numOfSegSites=sim.ALL_AVAIL, step=50),
        sim.PyEval(r"'Gen: %2d #seg sites: %d\n' % (gen, numOfSegSites)", step=50)
    ],
    gen = args.generations
)

logfile.write("Done simulating!\n")
logfile.write(time.strftime('%X %x %Z')+"\n")
logfile.write("----------\n")
logfile.flush()

locations = [pop.subPopIndPair(x)[0] for x in range(pop.popSize())]
rc.add_diploid_samples(nsamples=args.nsamples, 
                       sample_ids=pop.indInfo("ind_id"),
                       populations=locations)

del pop

logfile.write("Samples:\n")
logfile.write(str(rc.diploid_samples)+"\n")
logfile.write("----------\n")
logfile.flush()

rc.args.dump_sample_table(out=samples_file)

ts = rc.args.tree_sequence()
del rc

logfile.write("Loaded into tree sequence!\n")
logfile.write(time.strftime('%X %x %Z')+"\n")
logfile.write("----------\n")
logfile.flush()

# nodefile = open("nodes.txt","w")
# edgefile = open("edges.txt","w")
# ts.dump_text(nodes=nodefile, edgesets=edgefile)
# nodefile.flush()
# edgefile.flush()

minimal_ts = ts.simplify()
del ts

logfile.write("Simplified; now writing to treefile (if specified).\n")
logfile.write(time.strftime('%X %x %Z')+"\n")
logfile.write("----------\n")
logfile.flush()

if args.treefile is not None:
    minimal_ts.dump(args.treefile)

mut_seed=args.seed
logfile.write("Generating mutations with seed "+str(mut_seed)+"\n")
logfile.flush()

rng = msprime.RandomGenerator(mut_seed)
nodes = msprime.NodeTable()
edgesets = msprime.EdgesetTable()
sites = msprime.SiteTable()
mutations = msprime.MutationTable()
minimal_ts.dump_tables(nodes=nodes, edgesets=edgesets)
mutgen = msprime.MutationGenerator(rng, args.mut_rate)
mutgen.generate(nodes, edgesets, sites, mutations)

# print(nodes, file=logfile)
# print(edgesets, file=logfile)
# print(sites, file=logfile)
# print(mutations, file=logfile)

mutated_ts = msprime.load_tables(
    nodes=nodes, edgesets=edgesets, sites=sites, mutations=mutations)

del minimal_ts


logfile.write("Generated mutations!\n")
logfile.write(time.strftime('%X %x %Z')+"\n")
logfile.write("Mean pairwise diversity: {}\n".format(mutated_ts.get_pairwise_diversity()/mutated_ts.get_sequence_length()))
logfile.write("Sequence length: {}\n".format(mutated_ts.get_sequence_length()))
logfile.write("Number of trees: {}\n".format(mutated_ts.get_num_trees()))
logfile.write("Number of mutations: {}\n".format(mutated_ts.get_num_mutations()))

if args.outfile is None:
    print("NOT writing out genotype data.\n")
else:
    mutated_ts.write_vcf(outfile,ploidy=1)


logfile.write("All done!\n")
logfile.close()

