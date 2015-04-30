#!/usr/bin/python

##############
# GAUDIasm: Genetic Algorithms for Universal
# Design Inference and Atomic Scale Modeling
# Author: Jaime Rodriguez-Guerra Pedregal
# Email: jaime.rogue@gmail.com
# Web: https://bitbucket.org/jrgp/gaudi
#############

from gaudi.objectives import ObjectiveProvider
import random

def enable(**kwargs):
	return Dummy(**kwargs)

class Dummy(ObjectiveProvider):
	def __init__(self, *args, **kwargs):
		pass

	def evaluate(self, individual):
		return 1.0