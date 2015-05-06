#!/usr/bin/python

##############
# GAUDIasm: Genetic Algorithms for Universal
# Design Inference and Atomic Scale Modeling
# Author: Jaime Rodriguez-Guerra Pedregal
# Email: jaime.rogue@gmail.com
# Web: https://bitbucket.org/jrgp/gaudi
##############

from gaudi.genes import GeneProvider
from gaudi import box
import chimera
import deap.tools
import random

class Torsion(GeneProvider):
	"""
	Rotates bonds in given molecule. This plugins requires molecule.py.
	"""
	def __init__(self, target=None, flexibility=None, **kwargs):
		GeneProvider.__init__(self, **kwargs)
		self._kwargs = kwargs
		self.target = target
		self.flexibility = 360.0 if flexibility>360 else flexibility
		self.nonrotatable = ()
		try:
			self.compound = self.parent.genes[self.target].compound
		except KeyError:
			raise
		except AttributeError:
			raise
		else:
			self.rotatable_bonds = list(self.get_rotatable_bonds())
		
		self.allele = [self.random_angle() for b in self.rotatable_bonds]

	def express(self):
		for alpha, br in zip(self.allele, self.rotatable_bonds):
			try: 
				if all([a.idatmType in ('C2', 'N2') for a in br.bond.atoms]):
					alpha = 0 if alpha<180 else 180
				br.adjustAngle(alpha - br.angle, br.rotanchor)
			except AttributeError: # A null bondrot was returned -> non-rotatable bond
				continue
	
	def unexpress(self):
		for br in self.rotatable_bonds:
			br.reset()
	
	def mate(self, mate):
		self.allele[:], mate.allele[:] = deap.tools.cxSimulatedBinaryBounded(
				self.allele, mate.allele, eta=self.cxeta,
				low=-0.5*self.flexibility, up=0.5*self.flexibility)

	def mutate(self, indpb):
		self.allele, = deap.tools.mutPolynomialBounded(self.allele,
						indpb=self.indpb, eta=self.mteta, 
						low=-0.5*self.flexibility, 
						up=0.5*self.flexibility)

	def write(self, path, name):
		pass
	#####
	def random_angle(self):
		return random.uniform(-0.5*self.flexibility, 0.5*self.flexibility)

	def get_rotatable_bonds(self):
		bonds = set(b for a in self.compound.mol.atoms for b in a.bonds if not a.element.isMetal)
		bonds = sorted(bonds, key=lambda b: min(y.serialNumber for y in b.atoms))

		existing_bondrots_bonds = []
		for b in bonds:
			if b in existing_bondrots_bonds:
				continue
			a1,a2 = b.atoms
			if 	a1 not in self.nonrotatable and \
				a1.idatmType in ('C3', 'N3', 'C2', 'N2') and \
				(a1.numBonds > 1 and a2.numBonds > 1) or \
				a1.name == 'DUM' or a2.name == 'DUM' :
					try:
						br = chimera.BondRot(b)
					except (chimera.error, ValueError), v:
						if "cycle" in str(v):
							continue #discard bonds in cycles!
						elif "already used" in str(v):
							continue
						else:
							raise
					else:
						br.rotanchor = box.find_nearest(self.compound.donor, b.atoms)
						yield br

	def update_rotatable_bonds(self):
		self.rotatable_bonds = list(self.get_rotatable_bonds())


	def __deepcopy__(self, memo):
		new = self.__class__(self.target, self.flexibility, 
							**self._kwargs)
		new.__dict__.update((k,v) for k,v in self.__dict__.items())
		new.allele[:] = self.allele[:]
		return new 

def enable(**kwargs):
	return Torsion(**kwargs)