import msprime
import pyslim
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

for outdir in ["run_009673", "run_005464", "run_031486", "run_027034", "run_015598", "run_012720"]:

    treefile = outdir + "/results.trees"
    ts = pyslim.load(treefile, slim_format=True)

    def grid_samples(ts, n, m, prob=1.0):
        '''
        Returns a list of lists of node IDs, the ones that fall in each rectangle 
        given by discretizing into an (nxm) grid.
        '''
        samples = [[[] for _ in range(n)] for _ in range(m)]
        locs = np.array(ts.tables.individuals.location)
        locs.resize((int(len(ts.tables.individuals.location,)/3), 3))
        max_x = np.ceil(10*max(locs[:,0]))/10
        max_y = np.ceil(10*max(locs[:,1]))/10
        for ind in ts.individuals():
            if np.random.uniform() < prob:
                i = int(np.floor(ind.location[0] * n / max_x))
                j = int(np.floor(ind.location[1] * m / max_x))
                samples[i][j].extend(ind.nodes)
        return samples


    sample_grid = grid_samples(ts, 4, 4, 0.05)
    samples = [a for b in sample_grid for a in b]
    windows = np.linspace(0.0, ts.sequence_length, 1000)
    win_mids = (windows[1:] + windows[:-1])/2

    bs = msprime.BranchLengthStatCalculator(ts)

    divs = np.array(bs.divergence(samples, windows=windows))

    fig = plt.figure(figsize=(15,4))

    for i in range(divs.shape[1]):
        plt.plot(win_mids, divs[:,i])

    plt.savefig(outdir + "/divergences_1000.png", dpi=288)
