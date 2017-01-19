#!/usr/bin/env python3.5
description = '''
Convert simulations to msprime and VCF.
'''

import gzip
import sys
from optparse import OptionParser
import math
import time
import random
import msprime
import ftprime

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
parser.add_option("-i","--infile",dest="infile",help="input file (output by background-sim.py)")
parser.add_option("-k","--nsamples",dest="nsamples",help="number of *diploid* samples, total")
parser.add_option("-u","--mut_rate",dest="mut_rate",help="mutation rate",default=1e-7)
parser.add_option("-o","--outfile",dest="outfile",help="name of output file (or '-' for stdout)",default="-")
parser.add_option("-g","--logfile",dest="logfile",help="name of log file (or '-' for stdout)",default="-")
(options,args) =  parser.parse_args()

infile = fileopt(options.infile, "r")
outfile = fileopt(options.outfile, "w")
logfile = fileopt(options.logfile, "w")

logfile.write("Options:\n")
logfile.write(str(options)+"\n")

nsamples=int(options.nsamples)
mut_rate=float(options.mut_rate)

# Read in from the file::
#   generations
#   length
#   N

name,val = infile.readline().split(":")
if name != "# generations":
    raise ValueError("Bad file format.")
else:
    generations = int(val.strip())

name,val = infile.readline().split(":")
if name != "# length":
    raise ValueError("Bad file format.")
else:
    length = int(val.strip())

# total number of individuals; needed for translating indiv id to time
name,val = infile.readline().split(":")
if name != "# N":
    raise ValueError("Bad file format.")
else:
    N = int(val.strip())

def ind_to_time(k):
    return 1+generations-math.floor(k/N)

# Input is of this form:
# offspringID parentID startingPloidy rec1 rec2 ....
# ... coming in *pairs*

args=ftprime.ARGrecorder()
for k in range(N):
    for mapa in ftprime.mapa_labels:
        print("Adding"+str(ftprime.ind_to_chrom(k,mapa))+"\n")
        args.add_individual(ftprime.ind_to_chrom(k,mapa),ind_to_time(k))

while True:
    line=infile.readline()
    if not line:
        break
    print("A: "+line)
    child,parent,ploid,*rec = [int(x) for x in line.split()]
    for mapa in ftprime.mapa_labels:
        if mapa==ftprime.mapa_labels[1]:
            line=infile.readline()
            print(" : "+line)
            prev_child=child
            child,parent,ploid,*rec = [int(x) for x in line.split()]
            if child != prev_child:
                raise ValueError("Recombs not in pairs:"+str(child)+"=="+str(prev_child))
        start=0
        child_chrom=ftprime.ind_to_chrom(child,mapa)
        args.add_individual(child_chrom,ind_to_time(child))
        for r in rec:
            args.add_record(
                    left=start,
                    right=r,
                    parent=ftprime.ind_to_chrom(parent,ftprime.mapa_labels[ploid]),
                    children=(child_chrom,))
            start=r
            ploid=((ploid+1)%2)
        args.add_record(
                left=start,
                right=length,
                parent=ftprime.ind_to_chrom(parent,ftprime.mapa_labels[ploid]),
                children=(child_chrom,))


#######

# final pop IDs
pop_ids = range((generations-1)*N,generations*N)

samples=random.sample(pop_ids,nsamples)
# need chromosome ids
chrom_samples = [ ftprime.ind_to_chrom(x,a) for x in samples for a in ftprime.mapa_labels ]
meioser.records.add_samples(samples=chrom_samples)

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



