#!/usr/bin/env python3
description = '''
Simulate AND write to msprime/vcf.
'''

import gzip
import sys
from optparse import OptionParser
import math
import time
import random
from ftprime import RecombCollector, ind_to_chrom, mapa_labels
import msprime

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

parser = OptionParser(description=description)
parser.add_option("-T","--generations",dest="generations",help="number of generations to run for")
parser.add_option("-N","--popsize",dest="popsize",help="size of each subpopulation",default=100)
parser.add_option("-w","--gridwidth",dest="gridwidth",help="width of rectangular grid, in populations",default=3)
parser.add_option("-y","--gridheight",dest="gridheight",help="height of rectangular grid, in populations (default: equal to gridwidth)")
parser.add_option("-L","--length",dest="length",help="number of bp in the chromosome",default=1e4)
parser.add_option("-l","--nloci",dest="nloci",help="number of selected loci",default=20)
parser.add_option("-m","--migr",dest="migr",help="migration proportion between adjacent populations",default=.01)
parser.add_option("-u","--sel_mut_rate",dest="sel_mut_rate",help="mutation rate of selected alleles",default=1e-7)
parser.add_option("-r","--recomb_rate",dest="recomb_rate",help="recombination rate",default=2.5e-8)
parser.add_option("-S","--selection_coef",help="strength of selection",default=.1)
parser.add_option("-k","--nsamples",dest="nsamples",help="number of *diploid* samples, total")
parser.add_option("-A","--ancestor_age",dest="ancestor_age",help="time to ancestor above beginning of sim")
parser.add_option("-U","--mut_rate",dest="mut_rate",help="mutation rate",default=1e-7)
parser.add_option("-t","--treefile",dest="treefile",help="name of output file for trees (default: not output)",default=None)

parser.add_option("-o","--outfile",dest="outfile",help="name of output VCF file (default: not output)",default=None)
parser.add_option("-g","--logfile",dest="logfile",help="name of log file (or '-' for stdout)",default="-")
parser.add_option("-s","--selloci_file",dest="selloci_file",help="name of file to output selected locus information",default="sel_loci.txt")
parser.add_option("-e", "--samples_file", help="name of file to output information on samples (default=(dir)/samples.tsv)")

(options,args) =  parser.parse_args()

import simuOpt
simuOpt.setOptions(alleleType='mutant')
import simuPOP as sim
from simuPOP.demography import migr2DSteppingStoneRates, migrSteppingStoneRates

if options.outfile is not None:
    outfile = fileopt(options.outfile, "w")
logfile = fileopt(options.logfile, "w")
selloci_file = options.selloci_file
if args.samples_file is None:
    args.samples_file = os.path.join(os.path.dirname(args.treefile),"samples.tsv")
samples_file = fileopt(args.samples_file,"w")

logfile.write("Options:\n")
logfile.write(str(options)+"\n")
logfile.write(time.strftime('%X %x %Z')+"\n")
logfile.write("----------\n")
logfile.flush()

generations=int(options.generations)
popsize=int(options.popsize)
nloci=int(options.nloci)
gridwidth=int(options.gridwidth)
if options.gridheight is not None:
    gridheight=int(options.gridheight)
else:
    gridheight=gridwidth

selection_coef=float(options.selection_coef)
length=float(options.length)
recomb_rate=float(options.recomb_rate)
sel_mut_rate=float(options.sel_mut_rate)
mut_rate=float(options.mut_rate)
migr=float(options.migr)
nsamples=int(options.nsamples)
ancestor_age=float(options.ancestor_age)

npops=gridwidth*gridheight

# increase spacing between loci as we go along the chromosome
rel_positions=[0.0 for k in range(args.nloci-1)]
for k in range(1,args.nloci-1):
    rel_positions[k] = rel_positions[k-1] + random.expovariate(1)*(k**2)
pos_fac=length/(rel_positions[-1] + random.expovariate(1)*(k**2))
locus_position=[x*pos_fac for x in rel_positions] + [args.length]

# initially polymorphic alleles
init_freqs=[[k/100,1-k/100,0,0] for k in range(1,11)]
locus_classes=[min(len(init_freqs)-1,math.floor(random.expovariate(1))) for k in range(nloci)]
init_classes=[list(filter(lambda k: locus_classes[k]==x,range(nloci))) for x in range(len(init_freqs))]

logfile.write("Locus positions:\n")
logfile.write(str(locus_position)+"\n")
logfile.write("----------\n")
logfile.flush()


init_geno=[sim.InitGenotype(freq=init_freqs[k],loci=init_classes[k]) for k in range(len(init_freqs))]

# record recombinations
rc = RecombCollector(
        nsamples=nsamples, generations=generations, N=popsize*npops,
        ancestor_age=ancestor_age, length=length, locus_position=locus_position)

id_tagger = sim.IdTagger()
id_tagger.apply(pop)

# record recombinations
rc = RecombCollector(
        first_gen=pop.indInfo("ind_id"), ancestor_age=args.ancestor_age, 
                              length=args.length, locus_position=locus_position)
###
# modified from http://simupop.sourceforge.net/manual_svn/build/userGuide_ch5_sec9.html

class FixedFitness:
    def __init__(self,s):
        self.s = s
    def __call__(self, loc, alleles):
        # print(str(loc)+":"+str(alleles)+"\n")
        # needn't return fitness for alleles=(0,0) as simupop knows that's 1
        if 0 in alleles:
            return 1. - self.s
        else:
            return 1. - 2.*self.s


pop = sim.Population(
        size=[popsize]*npops, 
        loci=[nloci], 
        lociPos=locus_position,
        infoFields=['ind_id','fitness','migrate_to'])

if min(gridheight,gridwidth)==1:
    migr_rates=migrSteppingStoneRates(
        migr, n=max(gridwidth,gridheight), circular=False)
else:
    migr_rates=migr2DSteppingStoneRates(
        migr, m=gridwidth, n=gridheight, diagonal=False, circular=False)

pop.evolve(
    initOps=[
        sim.InitSex(),
    ]+init_geno,
    preOps=[
        sim.PyOperator(lambda pop: rc.increment_time() or True),
        sim.Migrator(
            rate=migr_rates,
            mode=sim.BY_PROBABILITY),
        sim.AcgtMutator(rate=[sel_mut_rate], model='JC69'),
        sim.PyMlSelector(FixedFitness(selection_coef),
            output=">>"+selloci_file),
    ],
    matingScheme=sim.RandomMating(
        ops=[
            id_tagger,
            sim.Recombinator(intensity=recomb_rate,
                output=rc.collect_recombs,
                infoFields="ind_id"),
        ] ),
    postOps=[
        sim.Stat(numOfSegSites=sim.ALL_AVAIL, step=50),
        sim.PyEval(r"'Gen: %2d #seg sites: %d\n' % (gen, numOfSegSites)", step=50)
    ],
    gen = generations
)


logfile.write("Done simulating!\n")
logfile.write(time.strftime('%X %x %Z')+"\n")
logfile.write("----------\n")
logfile.flush()

# writes out events in this form:
# offspringID parentID startingPloidy rec1 rec2 ....

locations = [pop.subPopIndPair(x)[0] for x in range(pop.popSize())]
rc.add_diploid_samples(nsamples=args.nsamples, sample_ids=pop.indInfo("ind_id"),
                       populations=locations)

del pop

logfile.write("Samples:\n")
logfile.write(str(rc.diploid_samples)+"\n")
logfile.write("----------\n")
logfile.flush()

rc.args.dump_sample_table(out=samples_file)

ts = rc.args.tree_sequence()

logfile.write("Loaded into tree sequence!\n")
logfile.write(time.strftime('%X %x %Z')+"\n")
logfile.write("----------\n")
logfile.flush()

minimal_ts = ts.simplify()
del ts

logfile.write("Simplified; now writing to treefile (if specified).\n")
logfile.write(time.strftime('%X %x %Z')+"\n")
logfile.write("----------\n")
logfile.flush()

mut_seed=random.randrange(1,1000)
logfile.write("Generating mutations with seed "+str(mut_seed)+"\n")
rng = msprime.RandomGenerator(mut_seed)
nodes = msprime.NodeTable()
edgesets = msprime.EdgesetTable()
sites = msprime.SiteTable()
mutations = msprime.MutationTable()
minimal_ts.dump_tables(nodes=nodes, edgesets=edgesets)
mutgen = msprime.MutationGenerator(rng, mut_rate)
mutgen.generate(nodes, edgesets, sites, mutations)
mutated_ts = msprime.load_tables(
    nodes=nodes, edgesets=edgesets, sites=sites, mutations=mutations)

del minimal_ts

if options.treefile is not None:
    mutated_ts.dump(options.treefile)

logfile.write("Generated mutations!\n")
logfile.write(time.strftime('%X %x %Z')+"\n")
logfile.write("Mean pairwise diversity: {}\n".format(mutated_ts.get_pairwise_diversity()/mutated_ts.get_sequence_length()))
logfile.write("Sequence length: {}\n".format(mutated_ts.get_sequence_length()))
logfile.write("Number of trees: {}\n".format(mutated_ts.get_num_trees()))
logfile.write("Number of mutations: {}\n".format(mutated_ts.get_num_mutations()))

if options.outfile is None:
    print("NOT writing out genotype data.\n")
else:
    mutated_ts.write_vcf(outfile,ploidy=1)


logfile.write("All done!\n")
logfile.close()
