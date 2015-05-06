#!/usr/bin/python

##############
# GAUDIasm: Genetic Algorithms for Universal
# Design Inference and Atomic Scale Modeling
# Author: Jaime Rodriguez-Guerra Pedregal
# Email: jaime.rogue@gmail.com
# Web: https://bitbucket.org/jrgp/gaudi
#############

from gaudi.objectives import ObjectiveProvider
import gaudi.parse
import chimera
import math

def enable(**kwargs):
	return Angle(**kwargs)

class Angle(ObjectiveProvider):
	def __init__(self, threshold=None, tolerance=-0.1, probes=None,
				*args, **kwargs):
		ObjectiveProvider.__init__(self, **kwargs)
		self.threshold = threshold
		self.tolerance = tolerance
		self.molecules = tuple(m.compound.mol for m in self.parent.genes 
							if m.__class__.__name__ == "Molecule")
		self._probes = probes
		self.probes = list(self._getprobes())
	
	def evaluate(self):
		atoms_coords = [a.molecule.openState.xform.apply(a.coord()) for a in self.probes]
		delta = 180.0
		try:
			angle = chimera.angle(*atoms_coords)
		except TypeError: # four atoms, means dihedral
			angle = chimera.dihedral(*atoms_coords)
		except TypeError: # threshold is str, calc abs sine
			if self.threshold == 'planar':
				delta = abs(math.sin(math.radians(angle)))
		else:
			delta = self.threshold - angle.real
		
		return delta

## TODO: Probes get lost if rotamers are applied!
	def _getprobes(self):
		for probe in self._probes:
			mol, serial = gaudi.parse.parse_rawstring(probe)
			try: 
				if isinstance(serial, int):
					atom = next(a for a in self.parent.genes[mol].compound.mol.atoms 
								if serial == a.serialNumber)
				else:
					atom = next(a for a in self.parent.genes[mol].compound.mol.atoms 
								if serial == a.name)
			except KeyError:
				print "Molecule not found"
				raise
			except StopIteration:
				print "No atoms matched for probe", probe
				raise
			else:
				yield atom