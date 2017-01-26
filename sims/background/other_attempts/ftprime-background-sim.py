#!/usr/bin/env python3.5
description = '''
Simulate.
'''

import gzip
import sys
from optparse import OptionParser
import math
import time
import msprime
import simuOpt
simuOpt.setOptions(alleleType='mutant')
import simuPOP as sim
import random
import ftprime
from itertools import accumulate

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
parser.add_option("-k","--nsamples",dest="nsamples",help="number of *diploid* samples, total")
parser.add_option("-N","--popsize",dest="popsize",help="size of each subpopulation",default=100)
parser.add_option("-w","--width",dest="width",help="width of square grid, in populations",default=3)
parser.add_option("-L","--length",dest="length",help="number of bp in the chromosome",default=1e4)
parser.add_option("-l","--nloci",dest="nloci",help="number of selected loci",default=20)
parser.add_option("-m","--migr",dest="migr",help="migration proportion between adjacent populations",default=.01)
parser.add_option("-u","--mut_rate",dest="mut_rate",help="mutation rate",default=1e-7)
parser.add_option("-r","--recomb_rate",dest="recomb_rate",help="recombination rate",default=2.5e-8)
parser.add_option("-a","--gamma_alpha",dest="gamma_alpha",help="alpha parameter in gamma distributed selection coefficients",default=1)
parser.add_option("-b","--gamma_beta",dest="gamma_beta",help="beta parameter in gamma distributed selection coefficients",default=.001)
parser.add_option("-o","--outfile",dest="outfile",help="name of output file (or '-' for stdout)",default="-")
parser.add_option("-g","--logfile",dest="logfile",help="name of log file (or '-' for stdout)",default="-")
(options,args) =  parser.parse_args()

outfile = fileopt(options.outfile, "w")
logfile = fileopt(options.logfile, "w")

logfile.write("Options:\n")
logfile.write(str(options)+"\n")

generations=int(options.generations)
nsamples=int(options.nsamples)
popsize=int(options.popsize)
nloci=int(options.nloci)
alpha=float(options.gamma_alpha)
beta=float(options.gamma_beta)
length=float(options.length)
recomb_rate=float(options.recomb_rate)
mut_rate=float(options.mut_rate)

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

###
# modified from http://simupop.sourceforge.net/manual_svn/build/userGuide_ch5_sec9.html

pop = sim.Population(
        size=[popsize]*50, 
        loci=[nloci], 
        lociPos=locus_position,
        infoFields=['ind_id','fitness'])
simu = sim.Simulator(pop)

meioser=ftprime.MeiosisTagger(nsamples,ngens=1+generations)
def step_gen(pop,param):
    # function to increment internal clock in MeiosisTagger
    # as in http://simupop.sourceforge.net/manual_svn/build/userGuide_ch5_sec14.html#python-operator-pyoperator
    dt,=param
    print("Time was", meioser.gen," and now is ",meioser.gen+dt)
    meioser.gen+=dt
    return True

reproduction = sim.RandomMating(
    ops=[
        sim.PyTagger(func=meioser.new_offspring),
        sim.Recombinator(intensity=recomb_rate),
    ]
)

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


simu.evolve(
    initOps=[
        sim.InitSex(),
        meioser
    ]+init_geno,
    preOps=[
        sim.AcgtMutator(rate=[0.0001], model='JC69'),
        sim.PyMlSelector(GammaDistributedFitness(alpha, beta),
            output='>>sel_loci.txt'),
        sim.PyOperator(func=step_gen,param=(1,)),
    ],
    matingScheme=reproduction,
    postOps=[
        sim.Stat(numOfSegSites=sim.ALL_AVAIL, step=50),
        sim.PyEval(r"'Gen: %2d #seg sites: %d\n' % (gen, numOfSegSites)", step=50)
    ],
    gen = generations
)

pop = simu.population(0)

pop_ids = [ ind.info('ind_id') for ind in pop.individuals() ]

samples=random.sample(pop_ids,nsamples)
# need chromosome ids
chrom_samples = [ ftprime.ind_to_chrom(x,a) for x in samples for a in ftprime.mapa_labels ]
times=[ meioser.time[x] for x in samples for a in ftprime.mapa_labels ]
meioser.records.add_samples(samples=chrom_samples,times=times,populations=[0 for x in chrom_samples])

print("msprime trees:")
ts=meioser.records.tree_sequence()
# ts.simplify()
ts.dump("sims.ts")

mut_seed=random.randrange(1,1000)
logfile.write("Generating mutations with seed "+str(mut_seed)+"\n")
rng = msprime.RandomGenerator(mut_seed)
ts.generate_mutations(mut_rate,rng)

logfile.write("Done!\n")
logfile.write(time.strftime('%X %x %Z')+"\n")
logfile.write("Mean pairwise diversity: {}\n".format(ts.get_pairwise_diversity()/ts.get_sequence_length()))
logfile.write("Sequence length: {}\n".format(ts.get_sequence_length()))
logfile.write("Number of trees: {}\n".format(ts.get_num_trees()))
logfile.write("Number of mutations: {}\n".format(ts.get_num_mutations()))

ts.write_vcf(outfile,ploidy=1)


