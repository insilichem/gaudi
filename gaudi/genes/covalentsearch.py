#!/usr/bin/python

##############
# GAUDIasm: Genetic Algorithms for Universal
# Design Inference and Atomic Scale Modeling
# Authors:  Jaime Rodriguez-Guerra Pedregal
#            <jaime.rodriguezguerra@uab.cat>
#           Jean-Didier Marechal
#            <jeandidier.marechal@uab.cat>
# Web: https://bitbucket.org/jrgp/gaudi
##############

"""
Document this!
"""

# GAUDI
from gaudi.genes import GeneProvider


class CovalentSearch(GeneProvider):

    def __init__(self):
        pass

    def express(self):
        pass

    def crossover(self, individual):
        pass

    def mutate(self):
        pass


def enable(**kwargs):
    return CovalentSearch(**kwargs)
