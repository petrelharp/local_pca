#!/usr/bin/python3.5

import sys, os
import time
import json
import random
import msprime

if len(sys.argv)>1:
    outdir=sys.argv[1]
else:
    print("Usage: {} (name of output directory)\n".format(sys.argv[0]))
    sys.exit(1)

if not os.path.isdir(outdir):
    os.mkdir(outdir)

base_options = {
        'Ne' : 1000,
        'm_rel' : 10,
        'M_rel' : 0.1,
        'T_rel' : 0.25,
        'nsamples' : 100,
        'chrom_len' : 25e6,
        'treefile' : "chrom{}.trees",
        'vcffile' : "chrom{}.vcf",
        'mut_rate' : 1e-6,
        'recomb_rate' : 1e-6,
     }

logfile=sys.stdout

options = { k:base_options.copy() for k in range(1,3) }
options[1]['Ne'] = base_options['Ne']/2
options[2]['Ne'] = base_options['Ne']*2

logfile.write("Begun: writing to {}\n".format(outdir))
logfile.write(time.strftime('%X %x %Z\n'))

with open(os.path.join(outdir,"config.json"),"w") as configfile:
    json.dump(options,configfile)

with open(os.path.join(outdir,"samples.tsv"),"w") as samplefile:
    samplefile.write("ID\tpopulation\n")
    for pop in range(1,5):
        for k in range(base_options['nsamples']):
            samplefile.write("{}\t{}\n".format(k,pop))

for chrom,opts in options.items():

    for fname in ('treefile','vcffile'):
        opts[fname]=os.path.join(outdir,opts[fname]).format(chrom)

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

    mut_seed=random.randrange(1,1000)
    rng = msprime.RandomGenerator(mut_seed)
    ts.generate_mutations(opts['mut_rate'],rng)

    ts.dump(opts['treefile'])

    ts.write_vcf(open(opts['vcffile'],'w'),ploidy=1)


logfile.write("Done!\n")
logfile.write(time.strftime('%X %x %Z\n'))
