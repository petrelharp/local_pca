#!/usr/bin/env python3
description = '''
Simulates on a rectangular grid that was bisected (in 'width') for some period in the past.
'''

import gzip
import sys, os
import math
import time
import random
from ftprime import RecombCollector
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
parser.add_argument("--post_generations","-R", type=int, dest="post_generations",
        help="number of generations since the division")
parser.add_argument("--split_generations","-S", type=int, dest="split_generations",
        help="number of generations of the split")
parser.add_argument("--pre_generations","-T", type=int, dest="pre_generations",
        help="number of generations to run before the split")
parser.add_argument("--popsize","-N", type=int, dest="popsize",
        help="size of each subpopulation",default=100)
parser.add_argument("--gridwidth","-w", type=int, dest="gridwidth",
        help="width of rectangular grid, in populations",default=3)
parser.add_argument("--gridheight","-y", type=int, dest="gridheight",
        help="height of rectangular grid, in populations (default: equal to gridwidth)")
parser.add_argument("--length","-L", type=float, dest="length",
        help="number of bp in the chromosome",default=1e4)
parser.add_argument("--nloci","-l", type=int, dest="nloci",
        help="number of selected loci",default=20)
parser.add_argument("--migr","-m", type=float, dest="migr",
        help="migration proportion between adjacent populations",default=.01)
parser.add_argument("--sel_mut_rate","-u", type=float, dest="sel_mut_rate",
        help="mutation rate of selected alleles",default=1e-7)
parser.add_argument("--recomb_rate","-r", type=float, dest="recomb_rate",
        help="recombination rate",default=2.5e-8)
parser.add_argument("--selection_coef","-s", type=float, dest="selection_coef",
	help="strength of selection",default=.1)
parser.add_argument("--nsamples","-k", type=int, dest="nsamples",
        help="number of *diploid* samples, total")
parser.add_argument("--ancestor_age","-A", type=float, dest="ancestor_age",
        help="time to ancestor above beginning of sim")
parser.add_argument("--mut_rate","-U", type=float, dest="mut_rate",
        help="mutation rate",default=1e-7)
parser.add_argument("--treefile","-t", type=str, dest="treefile",
        help="name of output file for trees (default: not output)",default=None)
parser.add_argument("--simplify_interval", "-I", dest="simplify_interval", type=int,
        help="Interval between simplify steps.", default=500)
parser.add_argument("--seed", "-d", dest="seed", type=int,
        help="random seed", default=random.randrange(1,1000))

parser.add_argument("--outfile","-o", type=str, dest="outfile",
        help="name of output VCF file (default: not output)",default=None)
parser.add_argument("--logfile","-g", type=str, dest="logfile",
        help="name of log file (or '-' for stdout)",default="-")
parser.add_argument("--selloci_file","-c", type=str, dest="selloci_file",
        help="name of file to output selected locus information",default="(dir)/sel_loci.txt")
parser.add_argument("--samples_file", "-e", type=str, dest="samples_file",
	help="name of file to output information on samples (default=(dir)/samples.tsv)")

args = parser.parse_args()

import simuOpt
import simuPOP as sim
# these just return the migration matrix as a list of lists
# from simuPOP.demography import migr2DSteppingStoneRates, migrSteppingStoneRates

def migrRates(migr, m, n, barrier=False):
    b = math.floor(m/2) - 1 # barrier is between rows b, b+1
    npops = m*n
    rows = [i for i in range(m) for j in range(n)]
    cols = [j for i in range(m) for j in range(n)]
    M = [[0.0 for u in range(npops)] for v in range(npops)]
    for u in range(npops):
        for v in range(npops):
            if ((abs(rows[u] - rows[v]) == 0 and abs(cols[u] - cols[v]) == 1)
                or (abs(rows[u] - rows[v]) == 1 and abs(cols[u] - cols[v]) == 0)):
                if not barrier or not (((rows[u] == b) and (rows[v] == b+1)) or ((rows[u] == b+1) and (rows[v] == b))):
                    M[u][v] = migr
        M[u][u] = 1.0 - sum(M[u])
    return M

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

if args.gridwidth < 2:
    raise ValueError("gridwidth must be at least 2 to bisect")

if args.gridheight is None:
	args.gridheight = args.gridwidth

logfile.write("Options:\n")
logfile.write(str(args)+"\n")
logfile.write(time.strftime('%X %x %Z')+"\n")
logfile.write("----------\n")
logfile.flush()

npops=args.gridwidth*args.gridheight

# loci evenly spaced on left half the choromosome
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

def fitness_fun(loc, alleles):
    if 0 in alleles:
        return 1. - args.selection_coef
    else:
        return 1. - 2. * args.selection_coef

pop = sim.Population(
        size=[args.popsize]*npops, 
        loci=[args.nloci,2], 
        lociPos=locus_position + [0, args.length],
        infoFields=['ind_id','fitness','migrate_to'])

id_tagger = sim.IdTagger()
id_tagger.apply(pop)

# set up recomb collector
first_gen = pop.indInfo("ind_id")
init_ts = msprime.simulate(2*len(first_gen),
                           length=max(locus_position))
haploid_labels = [(k,p) for k in first_gen 
                        for p in (0,1)]
node_ids = {x:j for x, j in zip(haploid_labels, init_ts.samples())}
rc = RecombCollector(ts=init_ts, node_ids=node_ids,
                     locus_position=locus_position)

migr_rates = migrRates(args.migr, m=args.gridwidth, n=args.gridheight, barrier=False)
barrier_rates = migrRates(args.migr, m=args.gridwidth, n=args.gridheight, barrier=True)

pop.evolve(
    initOps=[
        sim.InitSex(),
    ]+init_geno,
    preOps=[
        sim.PyOperator(lambda pop: rc.increment_time() or True),
        sim.Migrator(
            rate=migr_rates,
            mode=sim.BY_PROBABILITY,
            begin=0, 
            end=args.pre_generations),
        sim.Migrator(
            rate=barrier_rates,
            mode=sim.BY_PROBABILITY,
            begin=1 + args.pre_generations, 
            end=args.split_generations + args.pre_generations),
        sim.Migrator(
            rate=migr_rates,
            mode=sim.BY_PROBABILITY,
            begin=1 + args.split_generations + args.pre_generations,
            end=args.split_generations + args.pre_generations + args.post_generations),
        sim.SNPMutator(u=args.sel_mut_rate, v=args.sel_mut_rate),
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
        sim.PyEval(r"'Gen: %2d #seg sites: %d\n' % (gen, numOfSegSites)", step=50),
        sim.PyOperator(lambda pop: rc.simplify(pop.indInfo("ind_id")) or True, 
                       step=args.simplify_interval),
    ],
    gen = args.split_generations + args.pre_generations + args.post_generations
)

logfile.write("Done simulating!\n")
logfile.write(time.strftime('%X %x %Z')+"\n")
logfile.write("----------\n")
logfile.flush()

locations = [pop.subPopIndPair(x)[0] for x in range(pop.popSize())]
rc.add_locations(pop.indInfo("ind_id"), locations)

logfile.write("Collecting samples:\n")
logfile.write("  " + str(args.nsamples) + " of them")
logfile.write("  " + "ids:" + str(pop.indInfo("ind_id")))
logfile.write("  " + "locations:" + str(locations))

diploid_samples = random.sample(pop.indInfo("ind_id"), args.nsamples)
rc.simplify(diploid_samples)

del pop

logfile.write("Samples:\n")
logfile.write(str(rc.diploid_samples)+"\n")
logfile.write("----------\n")
logfile.flush()

ts = rc.args.tree_sequence()
del rc

logfile.write("Loaded into tree sequence!\n")
logfile.write(time.strftime('%X %x %Z')+"\n")
logfile.write("----------\n")
logfile.flush()

if args.treefile is not None:
    ts.dump(args.treefile)

logfile.write("Writing out samples.\n")
logfile.write(time.strftime('%X %x %Z')+"\n")
logfile.write("----------\n")
logfile.flush()

ts.dump_samples_text(samples_file)

mut_seed=args.seed
logfile.write(time.strftime('%X %x %Z')+"\n")
logfile.write("Generating mutations with seed "+str(mut_seed)+"\n")
logfile.flush()

rng = msprime.RandomGenerator(mut_seed)
nodes = msprime.NodeTable()
edgesets = msprime.EdgesetTable()
sites = msprime.SiteTable()
mutations = msprime.MutationTable()
ts.dump_tables(nodes=nodes, edgesets=edgesets)
mutgen = msprime.MutationGenerator(rng, args.mut_rate)
mutgen.generate(nodes, edgesets, sites, mutations)
mutated_ts = msprime.load_tables(
    nodes=nodes, edgesets=edgesets, sites=sites, mutations=mutations)

del ts

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
