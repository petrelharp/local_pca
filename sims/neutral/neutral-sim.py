#!/usr/bin/python3.5
description = '''
Simulate.
'''

import gzip
import sys, os
import math
import time
import random
import argparse

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

parser = argparse.ArgumentParser(description=description)
parser.add_argument("--nsamples", "-k", type=int, dest="nsamples", help="number of samples, total")
parser.add_argument("--popsize", "-N", type=int, dest="popsize", help="size of each subpopulation", default=100)
parser.add_argument("--width", "-w", type=int, dest="width", help="width of square grid, in populations", default=3)
parser.add_argument("--length", "-L", type=float, dest="length", help="number of bp in the chromosome", default=1e4)
parser.add_argument("--migr", "-m", type=float, dest="migr", help="migration proportion between adjacent populations", default=.01)
parser.add_argument("--min_recomb", "-r", type=float, dest="min_recomb", help="minimum recombination rate", default=2.5e-8)
parser.add_argument("--max_recomb", "-R", type=float, dest="max_recomb", help="maximum recombination rate", default=2.5e-8)
parser.add_argument("--mapfile", "-p", type=str, dest="mapfile", help="name of map file")
parser.add_argument("--mut_rate", "-u", type=float, dest="mut_rate", help="mutation rate", default=1e-7)
parser.add_argument("--tree_file", "-t", type=str, dest="tree_file", help="name of file to save tree sequence to.")
parser.add_argument("--samples_file", "-S", type=str, dest="samples_file", help="name of file to save sample information to.")
parser.add_argument("--outfile", "-o", type=str, dest="outfile", help="name of output file (or '-' for stdout).")
parser.add_argument("--logfile", "-g", type=str, dest="logfile", help="name of log file (or '-' for stdout)", default="-")
parser.add_argument("--seed", "-d", dest="seed", type=int,
        help="random seed", default=random.randrange(1,1000))

args = parser.parse_args()

if args.tree_file is None:
    print(description)
    raise(ValueError("Must specify output tree file."))

if args.outfile is None:
    args.outfile = os.path.join(os.path.dirname(args.tree_file),"sim.vcf")
if args.logfile is None:
    args.logfile = os.path.join(os.path.dirname(args.tree_file),"sim.log")
if args.samples_file is None:
    args.samples_file = os.path.join(os.path.dirname(args.tree_file),"samples.tsv")

logfile = fileopt(args.logfile, "w")

logfile.write("Options:\n")
logfile.write(str(args)+"\n")

outfile = fileopt(args.outfile, "w")
samples_file = fileopt(args.samples_file, "w")

sim.setRNG(seed=args.seed)
random.seed(args.seed)

popsize = float(args.popsize)
width = int(args.width)
nsamples = int(args.nsamples)
length = float(args.length)
migr = float(args.migr)
mut_rate = float(args.mut_rate)
if args.mapfile is not None:
    use_map = True
    mapfile = fileopt(args.mapfile, "r")
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
    min_recomb = float(args.min_recomb)
    max_recomb = float(args.max_recomb)
    use_map = (min_recomb!=max_recomb)
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

tree_sequence.dump_samples_text(samples_file)

if args.tree_file is not None:
    tree_sequence.dump(args.tree_file)

tree_sequence.write_vcf(outfile, ploidy=1)
logfile.close()

outfile.close()
logfile.close()

