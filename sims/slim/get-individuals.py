#!/usr/bin/env python3
description = '''
Get, and dump to text, the individuals table in the .trees output from a SLiM simulation.
'''

import sys, os
import gzip
import glob
import re
import argparse
import struct
import numpy as np

import msprime

parser = argparse.ArgumentParser(description=description)
parser.add_argument("--tree_file", "-t", type=str, nargs="*", dest="tree_file", 
                    help="name of file to load tree sequences from [default: .trees files in basedir.")
parser.add_argument("--basedir", "-o", type=str, dest="basedir", 
                    help="name of directory to save output files to.")
parser.add_argument("--indivfile", "-i", type=str, nargs="*", dest="indivfile", 
                    help="name of output files [default: as trees but with .indiv.tsv]")
parser.add_argument("--logfile", "-g", type=str, dest="logfile", 
                    help="name of log file")

args = parser.parse_args()
argdict = vars(args)

if args.basedir is None and args.indivfile is None:
    print(description)
    raise ValueError("Must specify at least basedir and indivfile (run with -h for help).")

if args.tree_file is None or len(args.tree_file) == 0:
    args.tree_file = glob.glob(os.path.join(args.basedir, "*.trees"))

if args.indivfile is None or len(args.indivfile) == 0:
    args.indivfile = [os.path.join(args.basedir, re.sub("[.]trees$", "", os.path.basename(x))) 
                        + ".indiv.tsv" for x in args.tree_file]

if args.logfile is None:
    args.logfile = os.path.join(args.basedir, "add_muts.log")

assert len(args.indivfile) == len(args.tree_file)

logfile = open(args.logfile, "w")

# typedef struct __attribute__((__packed__)) {
# 	slim_pedigreeid_t pedigree_id_;	// 8 bytes (int64_t): the SLiM pedigree ID for this individual, assigned by pedigree rec
# 	slim_age_t age_; // 4 bytes (int32_t)
# 	slim_objectid_t subpopulation_id_; // 4 bytes (int32_t)
# } IndividualMetadataRec;

class slimIndividual(object):
    def __init__(self, table_row):
        ped_id, age, subpop = struct.unpack("qii", table_row.metadata)
        location = table_row.location[:3]
        self.ped_id = ped_id
        self.age = age
        self.subpop = subpop
        self.location = location


header = ["id", "age", "subpop", "location_x", "location_y", "location_z", "node_ma", "node_pa"]

def write_indivs(chrom):
    treefile = args.tree_file[chrom]
    out = open(args.indivfile[chrom], "w")
    logfile.write("Reading trees from " + treefile + "\n")
    ts = msprime.load(treefile)
    node_inds = [np.where(ts.tables.nodes.individual == u)[0] 
                 for u in range(ts.tables.individuals.num_rows)]
    individuals = [slimIndividual(ts.tables.individuals[k]) 
                   for k in range(ts.tables.individuals.num_rows)]
    logfile.write("Saving to" + args.indivfile[chrom] + "\n")
    out.write("\t".join(header) + "\n");
    for k, ind in enumerate(individuals):
        data = [ind.ped_id, ind.age, ind.subpop] + list(ind.location) + list(node_inds[k])
        out.write("\t".join(map(str, data)) + "\n")
    out.close()
    return True

for k in range(len(args.tree_file)):
    write_indivs(k)

logfile.write("Done!\n")
logfile.close()

