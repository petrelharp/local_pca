#!/usr/bin/env python3.5
description = '''
Simulate AND write to msprime/vcf.
'''

import gzip
import sys
from optparse import OptionParser
import math
import time
import random
from itertools import accumulate
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
parser.add_option("-w","--width",dest="width",help="width of square grid, in populations",default=3)
parser.add_option("-L","--length",dest="length",help="number of bp in the chromosome",default=1e4)
parser.add_option("-l","--nloci",dest="nloci",help="number of selected loci",default=20)
parser.add_option("-m","--migr",dest="migr",help="migration proportion between adjacent populations",default=.01)
parser.add_option("-u","--sel_mut_rate",dest="sel_mut_rate",help="mutation rate of selected alleles",default=1e-7)
parser.add_option("-r","--recomb_rate",dest="recomb_rate",help="recombination rate",default=2.5e-8)
parser.add_option("-a","--gamma_alpha",dest="gamma_alpha",help="alpha parameter in gamma distributed selection coefficients",default=.23)
parser.add_option("-b","--gamma_beta",dest="gamma_beta",help="beta parameter in gamma distributed selection coefficients",default=5.34)
parser.add_option("-k","--nsamples",dest="nsamples",help="number of *diploid* samples, total")
parser.add_option("-A","--ancestor_age",dest="ancestor_age",help="time to ancestor above beginning of sim")
parser.add_option("-U","--mut_rate",dest="mut_rate",help="mutation rate",default=1e-7)
parser.add_option("-t","--treefile",dest="treefile",help="name of output file for trees (or '-' for stdout)")

parser.add_option("-o","--outfile",dest="outfile",help="name of output file (or '-' for stdout)",default="-")
parser.add_option("-g","--logfile",dest="logfile",help="name of log file (or '-' for stdout)",default="-")
parser.add_option("-s","--selloci_file",dest="selloci_file",help="name of file to output selected locus information",default="sel_loci.txt")
(options,args) =  parser.parse_args()

import simuOpt
simuOpt.setOptions(alleleType='mutant')
import simuPOP as sim
from simuPOP.demography import migr2DSteppingStoneRates

outfile = fileopt(options.outfile, "w")
logfile = fileopt(options.logfile, "w")
selloci_file = options.selloci_file

logfile.write("Options:\n")
logfile.write(str(options)+"\n")
logfile.write(time.strftime('%X %x %Z')+"\n")
logfile.write("----------\n")
logfile.flush()

generations=int(options.generations)
popsize=int(options.popsize)
nloci=int(options.nloci)
width=int(options.width)
alpha=float(options.gamma_alpha)
beta=float(options.gamma_beta)
length=float(options.length)
recomb_rate=float(options.recomb_rate)
sel_mut_rate=float(options.sel_mut_rate)
mut_rate=float(options.mut_rate)
migr=float(options.migr)
nsamples=int(options.nsamples)
ancestor_age=float(options.ancestor_age)

npops=width*width

# increase spacing between loci as we go along the chromosome
spacing_fac=9
rel_positions=list(accumulate([random.expovariate(1)*(1+spacing_fac*k/nloci) for k in range(nloci)]))
pos_fac=length/(rel_positions[-1]+random.expovariate(1)*(1+spacing_fac))
locus_position=[x*pos_fac for x in rel_positions]

# initially polymorphic alleles
init_freqs=[[k/100,1-k/100,0,0] for k in range(1,11)]
locus_classes=[min(len(init_freqs)-1,math.floor(random.expovariate(1))) for k in range(nloci)]
init_classes=[list(filter(lambda k: locus_classes[k]==x,range(nloci))) for x in range(len(init_freqs))]

init_geno=[sim.InitGenotype(freq=init_freqs[k],loci=init_classes[k]) for k in range(len(init_freqs))]

# record recombinations
rc = RecombCollector(
        nsamples=nsamples, generations=generations, N=popsize*npops,
        ancestor_age=ancestor_age, length=length, locus_position=locus_position)

###
# modified from http://simupop.sourceforge.net/manual_svn/build/userGuide_ch5_sec9.html

class GammaDistributedFitness:
    def __init__(self, alpha, beta):
        # mean is alpha/beta
        self.coefMap = {}
        self.alpha = alpha
        self.beta = beta
     
    def __call__(self, loc, alleles):
        # because s is assigned for each locus, we need to make sure the
        # same s is used for fitness of genotypes 01 (1-s) and 11 (1-2s)
        # at each locus
        if loc in self.coefMap:
            s = self.coefMap[loc]
        else:
            s = random.gammavariate(self.alpha, self.beta)
            self.coefMap[loc] = s
        # print(str(loc)+":"+str(alleles)+"\n")
        # needn't return fitness for alleles=(0,0) as simupop knows that's 1
        if 0 in alleles:
            return 1. - s
        else:
            return 1. - 2.*s

pop = sim.Population(
        size=[popsize]*npops, 
        loci=[nloci], 
        lociPos=locus_position,
        infoFields=['ind_id','fitness','migrate_to'])

pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.IdTagger(),
    ]+init_geno,
    preOps=[
        sim.Migrator(
            rate=migr2DSteppingStoneRates(
                migr, m=width, n=width, diagonal=False, circular=False),
            mode=sim.BY_PROBABILITY),
        sim.AcgtMutator(rate=[sel_mut_rate], model='JC69'),
        sim.PyMlSelector(GammaDistributedFitness(alpha, beta),
            output=">>"+selloci_file),
    ],
    matingScheme=sim.RandomMating(
        ops=[
            sim.IdTagger(),
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

rc.add_samples()

# not really the sample locs, but we don't care about that
sample_locs = [ (0,0) for _ in range(rc.nsamples) ]

ts = rc.args.tree_sequence(samples=sample_locs)

logfile.write("Loaded into tree sequence!\n")
logfile.write(time.strftime('%X %x %Z')+"\n")
logfile.write("----------\n")
logfile.flush()

ts.simplify(samples=list(range(nsamples)))

logfile.write("Simplified; now writing to treefile (if specified).\n")
logfile.write(time.strftime('%X %x %Z')+"\n")
logfile.write("----------\n")
logfile.flush()

mut_seed=random.randrange(1,1000)
logfile.write("Generating mutations with seed "+str(mut_seed)+"\n")
rng = msprime.RandomGenerator(mut_seed)
ts.generate_mutations(mut_rate,rng)

if options.treefile is not None:
    ts.dump(options.treefile)

logfile.write("Generated mutations!\n")
logfile.write(time.strftime('%X %x %Z')+"\n")
logfile.write("Mean pairwise diversity: {}\n".format(ts.get_pairwise_diversity()/ts.get_sequence_length()))
logfile.write("Sequence length: {}\n".format(ts.get_sequence_length()))
logfile.write("Number of trees: {}\n".format(ts.get_num_trees()))
logfile.write("Number of mutations: {}\n".format(ts.get_num_mutations()))

print("NOT writing out vcf to",outfile)
# ts.write_vcf(outfile,ploidy=1)


logfile.write("All done!\n")
logfile.close()
