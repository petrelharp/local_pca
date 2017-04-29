#!/usr/bin/env python3
description = '''
Simulate with simuPOP AND write to msprime/vcf:
    a population arrangement switches from A<-B|C to A|B->C 
    at some point in the past.

Output in selloci is of the form
   loc a1 a2 fitness
meaning locus number, allele1, allele2, fitness when it is first seen.
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
parser.add_argument('--relative_switch_time', '-w', default=0.25, type=float, 
        help="Time since rearrangement of populations in units of population size.")
parser.add_argument('--relative_fast_M', '-M', default=10, type=float, 
        help="Migration rate for 'close' pops in units of population size.")
parser.add_argument('--relative_slow_m', '-m', default=0.1, type=float, 
        help="Migration rate for 'distant' pops in units of population size.")
parser.add_argument("--generations", "-T", type=int, help="number of generations to run for before the switch")
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
simuOpt.setOptions(alleleType='mutant')
import simuPOP as sim
from simuPOP.demography import migr2DSteppingStoneRates, migrSteppingStoneRates

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
args.fast_M = args.relative_fast_M/(4*args.popsize)
args.slow_m = args.relative_slow_m/(4*args.popsize)
args.switch_time = int(args.relative_switch_time*(3*2*args.popsize))

logfile.write("Options:\n")
logfile.write(str(args)+"\n")
logfile.write(time.strftime('%X %x %Z')+"\n")
logfile.write("----------\n")
logfile.flush()

npops=3

# increase spacing between loci as we go along the chromosome
rel_positions=[0.0 for k in range(args.nloci)]
for k in range(1,args.nloci):
    rel_positions[k] = rel_positions[k-1] + random.expovariate(1)*(k**2)
pos_fac=args.length/(rel_positions[-1]+random.expovariate(1)*(args.nloci**2))
locus_position=[x*pos_fac for x in rel_positions]

# initially polymorphic alleles
init_freqs=[[k/100,1-k/100,0,0] for k in range(1,11)]
locus_classes=[min(len(init_freqs)-1,math.floor(random.expovariate(1))) for k in range(args.nloci)]
init_classes=[list(filter(lambda k: locus_classes[k]==x,range(args.nloci))) for x in range(len(init_freqs))]

logfile.write("Locus positions:\n")
logfile.write(str(locus_position)+"\n")
logfile.write("----------\n")
logfile.flush()

init_geno=[sim.InitGenotype(freq=init_freqs[k],loci=init_classes[k]) for k in range(len(init_freqs))]

def fitness_fun(loc, alleles):
    if 0 in alleles:
        return 1. - args.selection_coef
    else:
        return 1. - 2. * args.selection_coef

pop = sim.Population(
        size=[args.popsize]*npops, 
        loci=[args.nloci], 
        lociPos=locus_position,
        infoFields=['ind_id','fitness','migrate_to'])

id_tagger = sim.IdTagger()
id_tagger.apply(pop)

# record recombinations
rc = RecombCollector(
        first_gen=pop.indInfo("ind_id"), ancestor_age=args.ancestor_age, 
                              length=args.length, locus_position=locus_position)


migr_init = [ [ 0, args.slow_m, args.slow_m ],
              [ args.fast_M, 0, args.slow_m ],
              [ args.slow_m, args.slow_m, 0 ]]

migr_change = [ [ 0, args.slow_m, args.slow_m ],
                [ args.slow_m, 0, args.fast_M ],
                [ args.slow_m, args.slow_m, 0 ]]

# total number of generations to run simuPOP for
args.total_generations=args.generations+args.switch_time

pop.evolve(
    initOps=[
        sim.InitSex(),
    ]+init_geno,
    preOps=[
        sim.PyOperator(lambda pop: rc.increment_time() or True),
        sim.Migrator(
            rate=migr_init,
            mode=sim.BY_PROBABILITY,
            begin=0, end=args.switch_time),
        sim.Migrator(
            rate=migr_change,
            mode=sim.BY_PROBABILITY,
            begin=args.switch_time),
        sim.SNPmutator(u=args.sel_mut_rate, v=args.sel_mut_rate),
        sim.PyMlSelector(fitness_fun,
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
    gen = args.total_generations
)

logfile.write("Done simulating!\n")
logfile.write(time.strftime('%X %x %Z')+"\n")
logfile.write("----------\n")
logfile.flush()

locations = [pop.subPopIndPair(x)[0] for x in range(pop.popSize())]
rc.add_diploid_samples(pop.indInfo("ind_id"),locations)
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
