
# from http://simupop.sourceforge.net/manual_svn/build/userGuide_ch6_sec3.html
# to do our own recombination, create an OffspringGenerator,
# which gets called as a during-reproduction operator to create offspring,
# for instance:
def RandomMating(numOffspring=1., sexMode=RANDOM_SEX,
        ops=MendelianGenoTransmitter(), subPopSize=[],
        subPops=ALL_AVAIL, weight=0, selectionField='fitness'):
    'A basic diploid sexual random mating scheme.'
    return HomoMating(
        chooser=RandomParentsChooser(True, selectionField),
        generator=OffspringGenerator(ops, numOffspring, sexMode),
        subPopSize=subPopSize,
        subPops=subPops,
        weight=weight)

# In defining an OffspringGenerator...
# "Note that you need to specify all needed operators if you use parameter ops
# to change the operators used in a mating scheme (see Example
# HeteroMatingWeight). That is to say, you can use ops=Recombinator() to
# replace a default MendelianGenoTransmitter(), but you have to use
# ops=[IdTagger(), MendelianGenoTransmitter()] if you would like to add a
# during-mating operator to the default one."

# "A customized genotype transmitter is only a Python during-mating operator.
# Although it is possible to define a function and use a PyOperator directly
# (Example PyOperator), it is much better to derive an operator from
# PyOperator, as the case in Example newOperator:"

class sexSpecificRecombinator(PyOperator):
    def __init__(self, intensity=0, rates=0, loci=[], convMode=NO_CONVERSION,
            maleIntensity=0, maleRates=0, maleLoci=[], maleConvMode=NO_CONVERSION,
            *args, **kwargs):
        # This operator is used to recombine maternal chromosomes
        self.Recombinator = Recombinator(rates, intensity, loci, convMode)
        # This operator is used to recombine paternal chromosomes
        self.maleRecombinator = Recombinator(maleRates, maleIntensity,
            maleLoci, maleConvMode)
        #
        PyOperator.__init__(self, func=self.transmitGenotype, *args, **kwargs)
    #
    def transmitGenotype(self, pop, off, dad, mom):
        # Form the first homologous copy of offspring.
        self.Recombinator.transmitGenotype(mom, off, 0)
        # Form the second homologous copy of offspring.
        self.maleRecombinator.transmitGenotype(dad, off, 1)
        return True


# We'll need to use: http://simupop.sourceforge.net/manual_svn/build/refManual_ch3_sec5.html#class-genotransmitter
# class GenoTransmitter
#   This during mating operator is the base class of all genotype transmitters.
#   It is made available to users because it provides a few member functions
#   that can be used by derived transmitters, and by customized Python during
#   mating operators.
