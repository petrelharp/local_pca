#!/usr/bin/python3
usage = '''
Simulate a three-population neutral simulation with different Ne on different chromosomes.
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
parser.add_argument('--relative_fast_m', '-m', default=10, type=float, help="Migration rate for 'close' pops in units of Ne.")
parser.add_argument('--relative_slow_m', '-M', default=0.1, type=float, help="Migration rate for 'distant' pops in units of Ne.")

args = parser.parse_args()

if not os.path.isdir(args.outdir):
    os.mkdir(args.outdir)

base_options = {
        'Ne' : 1000,
        'm_rel' : args.relative_fast_m,
        'M_rel' : args.relative_slow_m,
        'T_rel' : args.relative_split_time,
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
    indiv=0
    for pop in range(1,5):
        for k in range(base_options['nsamples']):
            samplefile.write("msp_{}\tpop_{}\n".format(indiv,pop))
            indiv+=1

for chrom,opts in options.items():

    logfile.write("  simulating chromsome {}\n".format(chrom))
    logfile.write(time.strftime('     %X %x %Z\n'))
    logfile.flush()

    for fname in ('treefile','vcffile'):
        opts[fname]=os.path.join(args.outdir,opts[fname]).format(chrom)

    opts['m'] = opts['m_rel']/(4*opts['Ne'])
    opts['M'] = opts['M_rel']/(4*opts['Ne'])
    opts['T'] = opts['T_rel']*(4*opts['Ne'])

    npops = 3
    pops = [ msprime.PopulationConfiguration(sample_size=opts['nsamples'],
        initial_size=opts['Ne'], growth_rate=0.0)
        for _ in range(npops) ]

    migr_init = [ [ 0, opts['M'], opts['M'] ],
                  [ opts['m'], 0, opts['M'] ],
                  [ opts['M'], opts['M'], 0 ]]

    migr_change = [ msprime.MigrationRateChange(opts['T'], opts['m'], matrix_index=x) for x in [(1,2)] ] \
        + [ msprime.MigrationRateChange(opts['T'], opts['M'], matrix_index=x) for x in [(1,0)] ]

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

    rng = msprime.RandomGenerator(mut_seed)
    nodes = msprime.NodeTable()
    edgesets = msprime.EdgesetTable()
    sites = msprime.SiteTable()
    mutations = msprime.MutationTable()
    minimal_ts.dump_tables(nodes=nodes, edgesets=edgesets)
    mutgen = msprime.MutationGenerator(rng, opts['mut_rate'])
    mutgen.generate(nodes, edgesets, sites, mutations)

    mutated_ts = msprime.load_tables(
        nodes=nodes, edgesets=edgesets, sites=sites, mutations=mutations)

    del ts


    logfile.write("  done generating mutations! Writing out data.\n")
    logfile.write(time.strftime('     %X %x %Z\n'))
    logfile.flush()

    mutated_ts.dump(opts['treefile'])
    mutated_ts.write_vcf(open(opts['vcffile'],'w'),ploidy=1)


logfile.write("Done!\n")
logfile.write(time.strftime('     %X %x %Z\n'))
