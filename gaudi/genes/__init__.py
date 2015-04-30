#!/usr/bin/python

##############
# GAUDIasm: Genetic Algorithms for Universal
# Design Inference and Atomic Scale Modeling
# Author: Jaime Rodriguez-Guerra Pedregal
# Email: jaime.rogue@gmail.com
# Web: https://bitbucket.org/jrgp/gaudi
#############

import abc
from .. import plugin

class GeneProvider(object):
	"""
	From http://martyalchin.com/2008/jan/10/simple-plugin-framework/:
	Now that we have a mount point, we can start stacking plugins onto it. 
	As mentioned above, individual plugins will subclass the mount point. 
	Because that also means inheriting the metaclass, the act of subclassing 
	alone will suffice as plugin registration. Of course, the goal is to have
	plugins actually do something, so there would be more to it than just 
	defining a base class, but the point is that the entire contents of the
	class declaration can be specific to the plugin being written. The plugin 
	framework itself has absolutely no expectation for how you build the class, 
	allowing maximum flexibility. Duck typing at its finest.
	---
	So BaseGene is a mount point. Any gene implementing this should include these
	methods. Override them as needed! :)
	"""
	# This sole line is the magic behind the plugin system!
	__metaclass__ = plugin.PluginMount

	@abc.abstractmethod
	def express(self, individual):
		"""
		Compile the gene to an evaluable object.
		"""
	@abc.abstractmethod
	def unexpress(self, individual):
		"""
		Revert expression.
		"""

	@abc.abstractmethod
	def mutate(self):
		"""
		Perform a mutation on the gene.
		"""

	@abc.abstractmethod
	def crossover(self, gene):
		"""
		Perform a crossover with another gene of the same kind.
		"""

	@abc.abstractmethod
	def register(self):
		pass