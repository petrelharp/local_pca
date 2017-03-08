#!/usr/bin/python3.5
usage = '''
Simulate a bisplit neutral simulation with different Ne on different chromosomes.
'''

import sys, os
import time
import json
import random
import msprime

import argparse

parser = argparse.ArgumentParser(description=usage)
parser.add_argument('--outdir', '-o', help="Output directory.")
parser.add_argument('--Ne', '-N', nargs="+", type=int, help="List of effective population sizes used for each *chromosome*.")
parser.add_argument('--nsamples', '-n', type=int, help="Number of samples per population.")
parser.add_argument('--chrom_length', '-L', type=float, help="Length of chromosome in bp.")
parser.add_argument('--mut_rate', '-u', default=1e-7, type=float, help="Mutation rate per bp per generation.")
parser.add_argument('--recomb_rate', '-r', default=1e-7, type=float, help="Recombination rate per bp per generation.")
parser.add_argument('--relative_split_time', '-T', default=0.25, type=float, help="Time since rearrangement of populations in units of Ne.")

args = parser.parse_args()

if not os.path.isdir(args.outdir):
    os.mkdir(args.outdir)

base_options = {
        'Ne' : 1000,
        'm_rel' : 10,
        'M_rel' : 0.1,
        'T_rel' : 0.25,
        'nsamples' : args.nsamples,
        'chrom_len' : args.chrom_length,
        'treefile' : "chrom{}.trees",
        'vcffile' : "chrom{}.vcf",
        'mut_rate' : args.mut_rate,
        'recomb_rate' : args.recomb_rate,
     }

logfile=sys.stdout

options = { k:base_options.copy() for k in range(len(args.Ne)) }
for k in range(len(args.Ne)):
    options[k]['Ne'] = args.Ne[k]

logfile.write("Begun: writing to {}\n".format(args.outdir))
logfile.write(time.strftime('     %X %x %Z\n'))

with open(os.path.join(args.outdir,"config.json"),"w") as configfile:
    json.dump(options,configfile)

with open(os.path.join(args.outdir,"samples.tsv"),"w") as samplefile:
    samplefile.write("ID\tpopulation\n")
    for pop in range(1,5):
        for k in range(base_options['nsamples']):
            samplefile.write("{}\t{}\n".format(k,pop))

for chrom,opts in options.items():

    logfile.write("  simulating chromsome {}\n".format(chrom))
    logfile.write(time.strftime('     %X %x %Z\n'))
    logfile.flush()

    for fname in ('treefile','vcffile'):
        opts[fname]=os.path.join(args.outdir,opts[fname]).format(chrom)

    opts['m'] = opts['m_rel']/(4*opts['Ne'])
    opts['M'] = opts['M_rel']/(4*opts['Ne'])
    opts['T'] = opts['T_rel']*(4*opts['Ne'])

    pops = [ msprime.PopulationConfiguration(sample_size=opts['nsamples'],
        initial_size=opts['Ne'], growth_rate=0.0)
        for _ in range(4) ]

    migr_init = [ [ 0, opts['m'], opts['M'], opts['M'] ],
          [ opts['m'], 0, opts['M'], opts['M'] ],
          [ opts['M'], opts['M'], 0, opts['m'] ],
          [ opts['M'], opts['M'], opts['m'], 0 ]]

    migr_change = [ msprime.MigrationRateChange(opts['T'], opts['m'], matrix_index=x) for x in [(0,2),(1,3),(2,0),(3,1)] ] \
        + [ msprime.MigrationRateChange(opts['T'], opts['M'], matrix_index=x) for x in [(0,1),(0,3),(1,0),(1,2),(2,1),(2,3),(3,0),(3,2)] ]

    ts = msprime.simulate(
        Ne=opts['Ne'],
        length=opts['chrom_len'],
        recombination_rate=opts['recomb_rate'],
        population_configurations=pops,
        migration_matrix=migr_init,
        demographic_events=migr_change)

    logfile.write("  done simulating! Generating mutations.\n")
    logfile.write(time.strftime('     %X %x %Z\n'))
    logfile.flush()

    mut_seed=random.randrange(1,1000)
    rng = msprime.RandomGenerator(mut_seed)
    ts.generate_mutations(opts['mut_rate'],rng)

    logfile.write("  done generating mutations! Writing out data.\n")
    logfile.write(time.strftime('     %X %x %Z\n'))
    logfile.flush()

    ts.dump(opts['treefile'])
    ts.write_vcf(open(opts['vcffile'],'w'),ploidy=1)


logfile.write("Done!\n")
logfile.write(time.strftime('     %X %x %Z\n'))
