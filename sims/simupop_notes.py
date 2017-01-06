
###
# modified from http://simupop.sourceforge.net/manual_svn/build/userGuide_ch5_sec9.html

import simuOpt
simuOpt.setOptions(quiet=True, alleleType='mutant')
import simuPOP as sim
import random


pop = sim.Population(size=[100]*50, loci=[1000], infoFields=['fitness'])

class GammaDistributedFitness:
    def __init__(self, alpha, beta):
        self.coefMap = {}
        self.alpha = alpha
        self.beta = beta
     
    def __call__(self, loc, alleles):
        # because s is assigned for each locus, we need to make sure the
        # same s is used for fitness of genotypes 01 (1-s) and 11 (1-2s)
        # at each locus
        if loc in self.coefMap:
            s = self.coefMap[loc]
        else:
            s = random.gammavariate(self.alpha, self.beta)
            self.coefMap[loc] = s
        #
        if 0 in alleles:
            return 1. - s
        else:
            return 1. - 2.*s

pop.evolve(
    initOps=sim.InitSex(),
    preOps=[
        sim.AcgtMutator(rate=[0.00001], model='JC69'),
        sim.PyMlSelector(GammaDistributedFitness(0.23, 0.185),
            output='>>sel.txt'),
    ],
    matingScheme=sim.RandomMating(ops=sim.Recombinator(rates=0.1)),
    postOps=[
        sim.Stat(numOfSegSites=sim.ALL_AVAIL, step=50),
        sim.PyEval(r"'Gen: %2d #seg sites: %d\n' % (gen, numOfSegSites)",
            step=50)
    ],
    gen = 201
)
print(''.join(open('sel.txt').readlines()[:5]))
