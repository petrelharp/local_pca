#!/usr/bin/env python3
description = '''
Simulate a sequence of chromosomes, separately, using msprime, at different
parameter values.  The parameters that can vary between chromosomes are:
    - popsize
    - migration
    - min and max recombination rate
    - mutation rate
Each must be an argument either of length 1 or a list of lenth nchroms.

Samples will all be in the same configuration and order, with the same number
of samples from each subpopulation.
'''

import gzip
import sys, os
import math
import time
import random
import argparse
import multiprocessing

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
parser.add_argument("--nchroms", "-n", type=int, dest="nchroms", help="number of chromosomes")
parser.add_argument("--chrom_start", "-C", type=int, dest="chrom_start", help="index of first chromosome", default=0)
parser.add_argument("--nsamples", "-k", type=int, dest="nsamples", help="number of samples, total")
parser.add_argument("--popsize", "-N", type=float, nargs="*", dest="popsize", help="size of each subpopulation")
parser.add_argument("--width", "-w", type=int, dest="width", help="width of square grid, in populations")
parser.add_argument("--length", "-L", type=float, dest="length", help="number of bp in the chromosome", default=1e8)
parser.add_argument("--migr", "-m", type=float, nargs="*", dest="migr", help="migration proportion between adjacent populations")
parser.add_argument("--min_recomb", "-r", type=float, nargs="*", dest="min_recomb", help="minimum recombination rate", default=[2.5e-8])
parser.add_argument("--max_recomb", "-R", type=float, nargs="*", dest="max_recomb", help="maximum recombination rate", default=[2.5e-8])
parser.add_argument("--mapfile", "-p", type=str, dest="mapfile", help="name of map file")
parser.add_argument("--mut_rate", "-u", type=float, nargs="*", dest="mut_rate", help="mutation rate", default=[1e-8])
parser.add_argument("--basedir", "-o", type=str, dest="basedir", help="name of directory to save output files to. [default: msp_sim_$seed]")
parser.add_argument("--tree_file", "-t", type=str, dest="tree_file", help="name of file to save tree sequence to.")
parser.add_argument("--samples_file", "-S", type=str, dest="samples_file", help="name of file to save sample information to.")
parser.add_argument("--vcffile", "-v", type=str, dest="vcffile", help="name of VCF output file.")
parser.add_argument("--logfile", "-g", type=str, dest="logfile", help="name of log file")
parser.add_argument("--seed", "-d", dest="seed", type=int, help="random seed", default=random.randrange(1,1000))
parser.add_argument("--njobs", "-j", dest="njobs", type=int, help="number of parallel jobs", default=1)

args = parser.parse_args()

if args.basedir is None:
    args.basedir = "msp_sim_%04d" % args.seed

if not os.path.exists(args.basedir):
    os.makedirs(args.basedir)

if args.nsamples is None or args.nchroms is None or args.popsize is None or args.width is None or args.length is None:
    print(description)
    raise ValueError("Must specify more arguments (run with -h for help).")

argdict = vars(args)

vector_args = ['popsize', 'migr', 'min_recomb', 'max_recomb', 'mut_rate']
for a in vector_args:
    x = argdict[a]
    if x is not None and (len(x) == 1):
        argdict[a] = x * args.nchroms
    else: 
        if x is None or (len(x) != args.nchroms):
            raise ValueError(", ".join(vector_args) + "must all be of length 1 or of length nchroms.")

for a in ["tree_file", "logfile", "vcffile"]:
    x = argdict[a]
    if x is not None:
        if "%" not in x:
            raise ValueError("Output" + a + "must contain '%' to stand in for the chromosome number.")
        else:
            argdict[a].replace("%", "%02d")

if args.logfile is None:
    args.logfile = os.path.join(args.basedir, "sim.log")
if args.tree_file is None:
    args.tree_file = os.path.join(args.basedir, "sim%02d.trees")
if args.vcffile is None:
    args.vcffile = os.path.join(args.basedir, "sim%02d.vcf")
if args.samples_file is None:
    args.samples_file = os.path.join(args.basedir, "samples%02d.tsv")

logfile = fileopt(args.logfile, "w")

logfile.write("Options:\n")
logfile.write(str(args)+"\n")
logfile.flush()

random.seed(args.seed)
seeds = [random.randrange(1,1000) for _ in range(args.nchroms)]

def sim_chrom(chrom_num):

    random.seed(seeds[chrom_num])
    vcffile = fileopt(args.vcffile % chrom_num, "w")
    samples_file = fileopt(args.samples_file % chrom_num, "w")
    tree_file = args.tree_file % chrom_num

    popsize = float(args.popsize[chrom_num])
    width = int(args.width)
    nsamples = int(args.nsamples)
    length = float(args.length)
    migr = float(args.migr[chrom_num])
    mut_rate = float(args.mut_rate[chrom_num])
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
        min_recomb = float(args.min_recomb[chrom_num])
        max_recomb = float(args.max_recomb[chrom_num])
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

    msp_args = {'Ne' : popsize, 'mutation_rate' : mut_rate, 'length' : length}
    if width == 1:
        msp_args['sample_size'] = nsamples
    else:
        msp_args['population_configurations'] = pop_config
        msp_args['migration_matrix'] = mig_mat
    if use_map:
        msp_args['recombination_map'] = recomb_map
    else:
        msp_args['recombination_rate'] = max_recomb

    tree_sequence = msprime.simulate(**msp_args)

    logfile.write("Done!\n")
    logfile.write(time.strftime('%X %x %Z')+"\n")
    logfile.write("Mean pairwise diversity: {}\n".format(tree_sequence.get_pairwise_diversity()/length))
    logfile.write("Number of trees: {}\n".format(tree_sequence.get_num_trees()))
    logfile.write("Number of mutations: {}\n".format(tree_sequence.get_num_mutations()))
    logfile.flush()

    tree_sequence.dump_samples_text(samples_file)

    tree_sequence.dump(tree_file)

    tree_sequence.write_vcf(vcffile, ploidy=1)

    samples_file.close()
    vcffile.close()

    return True

if args.njobs > 1:
    p = multiprocessing.Pool(args.njobs)
    p.map(sim_chrom, range(args.chrom_start, args.chrom_start+args.nchroms))
else:
    for j in range(args.chrom_start, args.chrom_start+args.nchroms):
        sim_chrom(j)

logfile.close()

