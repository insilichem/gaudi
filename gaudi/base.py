#!/usr/bin/python

##############
# GAUDIasm: Genetic Algorithms for Universal
# Design Inference and Atomic Scale Modeling
# Author: Jaime Rodriguez-Guerra Pedregal
# Email: <jaime.rogue@gmail.com>
# Web: https://bitbucket.org/jrgp/gaudi
##############

from collections import OrderedDict
from zipfile import ZipFile, ZIP_DEFLATED, ZIP_STORED
import os
import pprint

import deap.base
import yaml

import gaudi.plugin

pp = pprint.PrettyPrinter(4)
class Individual(object):
	""" 
	Each individual is a potential solution. 
	It contains all that is needed for an evaluation. With multiprocessing in mind,
	individuals should be self-contained so it can be passed between threads.
	This replaces the need for all that initialization scripting in launch.py.

	:genes:    A dict containing parsed features (genes) for the solution
	:environment:   Constants of the system
	"""
	_CACHE = {}
	def __init__(self, genescfg=None, cache=None, cfg=None):
		self._genescfg = genescfg
		self.genes = OrderedDict()
		self.cfg = cfg
		gaudi.plugin.load_plugins(genescfg, container=self.genes, 
								cache=self._CACHE, parent=self)

	def evaluate(self):
		self.express()
		score = self.fitness.evaluate(self)
		self.unexpress()
		return score

	def express(self):
		""" 
		Express genes in this environment. Very much like 'compiling' the
		individual to a chimera.Molecule.
		"""
		for name, gene in self.genes.items():
			print "Expressing gene", name, "with allele"
			pp.pprint(gene.allele)
			gene.express()

	def unexpress(self):
		""" 
		Undo .express()
		"""
		for gene in reversed(self.genes.values()):
			gene.unexpress()

	def mate(self, individual):
		for gene in self.genes.values():
			gene.mate(individual.genes[gene.name])

		return self, individual

	def mutate(self, indpb):
		for gene in self.genes.values():
			gene.mutate(indpb)
		return self,

	def similar(self, individual):
		pass
	
	def write(self, path, name, i, compress=True):
		"""
		# Maybe someday we can pickle it all :/
		filename = os.path.join(path, '{}_{}.pickle.gz'.format(name,i))
		with gzip.GzipFile(filename, 'wb') as f:
			cPickle.dump(self, f, 0)
		return filename
		"""
		COMPRESS = ZIP_DEFLATED if compress else ZIP_STORED
		self.express()
		zipfilename = os.path.join(path, '{}__{:03d}.zip'.format(name,i))
		with ZipFile(zipfilename, 'w', COMPRESS) as z:
			output = OrderedDict() 
			for gene in self.genes.values():
				print "Writing", gene.name
				filename = gene.write(path, name)
				if filename:
					z.write(filename, os.path.basename(filename))
					os.remove(filename)
					output[gene.name] = os.path.basename(filename)
			try:
				output['score'] = list(self.fitness.values)
			except AttributeError: #fitness not in individual :/
				raise
			z.writestr('{}__{:03d}.gaudi'.format(name,i),
						yaml.dump(output, default_flow_style=False))
		self.unexpress()
		return zipfilename

class Fitness(deap.base.Fitness):
	"""
	Subclass of deap's Fitness to include details of objectives being evaluated
	and a function to evaluate them all at once. Since Fitness it's an Attribute
	of every individual, it should result in a self-contained object.
	"""
	objectives = OrderedDict()
	def __init__(self, parent=None, *args, **kwargs):
		deap.base.Fitness.__init__(self, *args, **kwargs)
		self.parent = parent
		if not self.objectives:
			gaudi.plugin.load_plugins(self.objectivelist, container=self.objectives)

	def evaluate(self, individual):
		score = []
		for name,obj in self.objectives.items():
			print "Evaluating", name
			score.append(obj.evaluate(individual))
		return score