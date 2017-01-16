#!/usr/bin/python3.5
description = '''
Simulate.
'''

import gzip
import sys
from optparse import OptionParser
import math
import time

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
parser.add_option("-k","--nsamples",dest="nsamples",help="number of samples, total")
parser.add_option("-N","--popsize",dest="popsize",help="size of each subpopulation",default=100)
parser.add_option("-w","--width",dest="width",help="width of square grid, in populations",default=3)
parser.add_option("-L","--length",dest="length",help="number of bp in the chromosome",default=1e4)
parser.add_option("-m","--migr",dest="migr",help="migration proportion between adjacent populations",default=.01)
parser.add_option("-r","--min_recomb",dest="min_recomb",help="minimum recombination rate",default=2.5e-8)
parser.add_option("-R","--max_recomb",dest="max_recomb",help="maximum recombination rate",default=2.5e-8)
parser.add_option("-p","--mapfile",dest="mapfile",help="name of map file")
parser.add_option("-u","--mut_rate",dest="mut_rate",help="mutation rate",default=1e-7)
parser.add_option("-o","--outfile",dest="outfile",help="name of output file (or '-' for stdout)",default="-")
parser.add_option("-g","--logfile",dest="logfile",help="name of log file (or '-' for stdout)",default="-")
(options,args) =  parser.parse_args()

outfile = fileopt(options.outfile, "w")
logfile = fileopt(options.logfile, "w")

logfile.write("Options:\n")
logfile.write(str(options)+"\n")

popsize = float(options.popsize)
width = int(options.width)
nsamples = int(options.nsamples)
length = float(options.length)
migr = float(options.migr)
mut_rate = float(options.mut_rate)
if options.mapfile is not None:
    use_map = True
    mapfile = fileopt(options.mapfile, "r")
    positions=[]
    rates=[]
    mapfile.readline()
    for line in mapfile:
        pos, rate, = map(float, line.split()[1:3])
        positions.append(pos)
        # Rate is expressed in centimorgans per megabase, which
        # we convert to per-base rates
        rates.append(rate * 1e-8)
    recomb_map = msprime.RecombinationMap(positions,rates)
else:
    min_recomb = float(options.min_recomb)
    max_recomb = float(options.max_recomb)
    use_map = (min_recomb==max_recomb)
    # recombination map
    n_recomb_steps = 100
    recomb_map = msprime.RecombinationMap(
            positions=[ k*length/n_recomb_steps for k in range(n_recomb_steps+1) ],
            rates=[ max_recomb - (max_recomb-min_recomb)*k/n_recomb_steps for k in range(n_recomb_steps+1) ])

logfile.write(time.strftime('%X %x %Z')+"\n")
logfile.write("Beginning simulation:\n")
logfile.flush()

# population setup
per_samples = math.ceil(nsamples/width/width)
pop_config = [ msprime.PopulationConfiguration(sample_size=per_samples) for i in range(width) for j in range(width) ]
mig_mat = [ [ migr if abs(i-k)+abs(j-l)==1 else 0.0 for i in range(width) for j in range(width) ] for k in range(width) for l in range(width) ]

logfile.write("Sample configuration:\n")
logfile.write(str([x.sample_size for x in pop_config])+"\n")
# logfile.write("Migration matrix:\n")
# logfile.write(str(mig_mat)+"\n")

if width==1:
    if not use_map:
        tree_sequence = msprime.simulate(
                sample_size=nsamples, 
                Ne=popsize,
                length=length, 
                recombination_rate=max_recomb,
                mutation_rate=mut_rate)
    else:
        tree_sequence = msprime.simulate(
                sample_size=nsamples, 
                Ne=popsize,
                recombination_map=recomb_map,
                mutation_rate=mut_rate)

else:
    if not use_map:
        tree_sequence = msprime.simulate(
                population_configurations=pop_config,
                migration_matrix=mig_mat,
                Ne=popsize,
                length=length,
                recombination_rate=max_recomb,
                mutation_rate=mut_rate)
    else:
        tree_sequence = msprime.simulate(
                population_configurations=pop_config,
                migration_matrix=mig_mat,
                Ne=popsize,
                recombination_map=recomb_map,
                mutation_rate=mut_rate)

logfile.write("Done!\n")
logfile.write(time.strftime('%X %x %Z')+"\n")
logfile.write("Mean pairwise diversity: {}\n".format(tree_sequence.get_pairwise_diversity()/length))
logfile.write("Number of trees: {}\n".format(tree_sequence.get_num_trees()))
logfile.write("Number of mutations: {}\n".format(tree_sequence.get_num_mutations()))
logfile.flush()

tree_sequence.write_vcf(outfile,ploidy=1)
logfile.close()

outfile.close()
logfile.close()

