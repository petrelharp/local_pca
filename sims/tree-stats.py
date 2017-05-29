#!/usr/bin/env python3
description = '''
Compute mean divergence between populations in windows along the genome.
'''

import gzip, csv
import sys, os, math
import msprime
import argparse

parser = argparse.ArgumentParser(description=description)
parser.add_argument('--treefile', '-t', type=str, 
        help="Name of the tree sequence file output by msprime.")
parser.add_argument('--window_size', '-w', type=float, 
        help="Width of the (regularly spaced) windows.")
parser.add_argument('--n_window', '-n', type=int, 
        help="Number of (regularly spaced) windows.")
parser.add_argument('--outfile', '-o', type=str, 
        help="Output file name.")

args = parser.parse_args()

ts = msprime.load(args.treefile)

pops = {}
ids = {} # record msprime -> ftprime ID here (not used)

for node_id in ts.samples():
    node = ts.node(node_id)
    pop = str(node.population)
    if not pop in pops:
        pops[pop] = []
    pops[pop] += [node_id]

leaf_sets = [list(map(int,u)) for u in pops.values()]
pop_names = list(pops.keys())

print('leaf sets:', leaf_sets)
print('samples:', list(ts.samples()))
print('ids:', list(ids.values()))

output_names = ["_".join([pop_names[i],pop_names[j]]) for i in range(len(pop_names))
                        for j in range(i,len(pop_names))]

if args.n_window is not None:
    if args.window_size is not None:
        raise ValueError("Can't specify both number of windows and window size.")
    else:
        args.window_size = ts.sequence_length / args.n_window
else:
    args.n_window = math.ceil(ts.sequence_length / args.window_size)

windows = [args.window_size * k for k in range(args.n_window + 1)]

tmrcas = ts.get_mean_tmrca(leaf_sets, windows)

with open(args.outfile, "w", newline="") as outfile:
    writer = csv.writer(outfile, delimiter="\t")
    writer.writerow(["start", "end"] + output_names)
    for k in range(len(tmrcas)):
        writer.writerow([windows[k], windows[k+1]] + tmrcas[k])

outfile.close()
