#!/usr/bin/python

# gaudi
# Genetic Algorithm for Unified Docking Inference
# A docking module for UCSF Chimera
# Jaime RGP <https://bitbucket.org/jrgp> @ UAB, 2014
import os
import yaml

class Settings(object):
	''' Simple parser for YAML settings file. '''
	def __init__(self, path, asDict=False):
		self._path = path
		self.data = {}
		self._parse()
	
	def weights(self):
		w = []
		for obj in self.objectives:
			if obj.type == 'hbonds' and hasattr(obj, 'targets') and len(obj.targets):
				w.append(obj.weight)
			w.append(obj.weight)
		return w
	
	def list_objectives(self):
		o = []
		for obj in self.objectives:
			o.append(obj.name)
			if obj.type == 'hbonds' and hasattr(obj, 'targets') and len(obj.targets):
				o.append('pref_'+obj.name)	
		return o

	def _parse(self):
		with open(self._path, 'r') as f:
			self.data = yaml.load(f)
		for k,v in self.data.items():
			if isinstance(v, list):
				self.__dict__[k] = [Param(d) for d in v]
			else:
				self.__dict__[k] = Param(v)

class Param(object):
	def __init__(self, *d):
		for d_ in d:
			self.__dict__.update(d_)

	def __str__(self):
		return "{}: {}-type objective with weight={}".format(self.name, self.type, self.weight)

#####
def _test_rebuild(cfg):
	for s, c in cfg.items():
		if isinstance(c, dict):
			print '\n['+s+']'
			for k, v in c.items():
				print k, "=", v
		elif isinstance(c, list):
			for i, l in enumerate(c):
				print '\n['+str(s)+' '+str(i)+']'
				for k, v in l.items():
					print k, "=", v

if __name__ == '__main__':
	import sys
	cfg = Settings(sys.argv[1])
	print [o.type for o in cfg.objective]

