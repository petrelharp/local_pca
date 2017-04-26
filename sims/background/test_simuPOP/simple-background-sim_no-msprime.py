#!/usr/bin/env python3
description = '''
Simulate a single population with varying levels of background selection.
'''

import gzip
import sys, os
import math
import time
import random

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
parser.add_argument("--generations","-T", type=int, help="number of generations to run for")
parser.add_argument("--popsize","-N", type=int, help="size of each subpopulation",default=100)
parser.add_argument("--length","-L", type=float, help="number of bp in the chromosome",default=1e4)
parser.add_argument("--nloci","-l", type=int, help="number of selected loci",default=20)
parser.add_argument("--sel_mut_rate","-u", type=float, help="mutation rate of selected alleles",default=1e-7)
parser.add_argument("--recomb_rate","-r", type=float, help="recombination rate",default=2.5e-8)
parser.add_argument("--selection_coef","-S", type=float, help="strength of selection",default=.1)
parser.add_argument("--nsamples","-k", type=int, help="number of *diploid* samples, total")
parser.add_argument("--ancestor_age","-A", type=int, help="time to ancestor above beginning of sim")
parser.add_argument("--mut_rate","-U", type=float, help="mutation rate",default=1e-7)
parser.add_argument("--treefile","-t", type=str, help="name of output file for trees (default: not output)")

parser.add_argument("--outfile","-o", type=str, help="name of output VCF file (default: not output)")
parser.add_argument("--logfile","-g", type=str, help="name of log file (or '-' for stdout)",default="-")
parser.add_argument("--selloci_file","-s", type=str, help="name of file to output selected locus information")
parser.add_argument("--samples_file", "-e", help="name of file to output information on samples (default=(dir)/samples.tsv)")

args = parser.parse_args()

import simuOpt
simuOpt.setOptions(alleleType='mutant')
import simuPOP as sim

if args.selloci_file is None:
    args.selloci_file = os.path.join(os.path.dirname(args.treefile),"sel_loci.txt")
if args.samples_file is None:
    args.samples_file = os.path.join(os.path.dirname(args.treefile),"samples.tsv")
samples_file = fileopt(args.samples_file,"w")

if args.outfile is not None:
    outfile = fileopt(args.outfile, "w")
logfile = fileopt(args.logfile, "w")

logfile.write("Options:\n")
logfile.write(str(args)+"\n")
logfile.write(time.strftime('%X %x %Z')+"\n")
logfile.write("----------\n")
logfile.flush()

# increase spacing between loci as we go along the chromosome
rel_positions=[0.0 for k in range(args.nloci)]
for k in range(args.nloci):
    rel_positions[k] = rel_positions[k-1] + random.expovariate(1)
pos_fac=args.length/(rel_positions[-1]+random.expovariate(1))
locus_position=[x*pos_fac for x in rel_positions]

# initially polymorphic alleles
init_freqs=[[k/100,1-k/100,0,0] for k in range(10)]
locus_classes=[min(len(init_freqs)-1,math.floor(random.expovariate(1))) for k in range(args.nloci)]
init_classes=[list(filter(lambda k: locus_classes[k]==x,range(args.nloci))) for x in range(len(init_freqs))]

logfile.write("Locus positions:\n")
logfile.write(str(locus_position)+"\n")
logfile.write("----------\n")
logfile.flush()


init_geno=[sim.InitGenotype(freq=init_freqs[k],loci=init_classes[k]) for k in range(len(init_freqs))]

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
        size=args.popsize,
        loci=[args.nloci], 
        lociPos=locus_position,
        infoFields=['ind_id','fitness','migrate_to'])


pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.IdTagger(),
    ]+init_geno,
    preOps=[
        sim.SNPMutator(u=args.sel_mut_rate, v=args.sel_mut_rate),
        sim.PyMlSelector(FixedFitness(args.selection_coef),
            output=">>"+args.selloci_file),
    ],
    matingScheme=sim.RandomMating(
        ops=[
            sim.IdTagger(),
            sim.Recombinator(intensity=args.recomb_rate,
                infoFields="ind_id"),
        ] ),
    postOps=[
        sim.Stat(numOfSegSites=sim.ALL_AVAIL, step=10),
        sim.PyEval(r"'Gen: %2d #seg sites: %d\n' % (gen, numOfSegSites)", step=50)
    ],
    gen = args.generations
)

logfile.write("Done simulating!\n")
logfile.write(time.strftime('%X %x %Z')+"\n")
logfile.write("----------\n")
logfile.flush()


logfile.write("All done!\n")
logfile.close()
