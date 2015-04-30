#!/usr/bin/python

##############
# GAUDIasm: Genetic Algorithms for Universal
# Design Inference and Atomic Scale Modeling
# Author: Jaime Rodriguez-Guerra Pedregal
# Email: jaime.rogue@gmail.com
# Web: https://bitbucket.org/jrgp/gaudi
##############

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