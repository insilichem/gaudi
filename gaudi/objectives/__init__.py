#!/usr/bin/python

##############
# GAUDIasm: Genetic Algorithms for Universal
# Design Inference and Atomic Scale Modeling
# Author: Jaime Rodriguez-Guerra Pedregal
# Email: jaime.rogue@gmail.com
# Web: https://bitbucket.org/jrgp/gaudi
##############
import chem, geom
import abc
from .. import plugin


class ObjectiveProvider(object):
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
	Mount point for plugins implementing new objectives to be evaluated by DEAP.
 	An objective only needs to accept a gaudi.Individual and return its score
 	with an evaluate() function. Apart from that, there's no requirements.
 	"""
	__metaclass__ = plugin.PluginMount
	def __init__(self, parent=None, name=None, weight=None, cache=None, environment=None,
				**kwargs):
		self.parent = parent
		self.name = name
		self.weight = weight
		self.env = environment
		try:
			self._cache = cache[self.name]
		except KeyError:
			self._cache = cache[self.name] = {}

	@abc.abstractmethod
	def evaluate(self, individual):
		"""
		Return the score of the individual under the current conditions.
		"""

