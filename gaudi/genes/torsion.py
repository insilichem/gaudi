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
	def __init__(self, parent=None, target=None, flexibility=None, **kwargs):
		self.parent = parent
		self.target = target
		self.flexibility = 360.0 if flexibility>360 else flexibility
		try:
			self.molecule = self.parent.genes[self.target]
		except KeyError:
			raise
		else:
			self.rotatable_bonds = self.get_rotatable_bonds()
		finally:
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
				self.allele, mate.allele,
				low=-0.5*self.flexibility, up=0.5*self.flexibility)

	def mutate(self):
		self.allele[:] = deap.tools.mutPolynomialBounded(self.allele, 
						low=-0.5*self.flexibility, 
						up=0.5*self.flexibility)[0]

	#####
	def random_angle(self):
		return random.uniform(-0.5*self.flexibility, 0.5*self.flexibility)

	def get_rotatable_bonds(self):
		try:
			existing_bondrots = self.rotatable_bonds
		except AttributeError:
			self.rotatable_bonds = []
		else:
			existing_bondrots_bonds = []
			for br in existing_bondrots:
				existing_bondrots_bonds.append(br.bond)
				yield br

		bonds = set(b for a in self.mol.atoms for b in a.bonds if not a.element.isMetal)
		bonds = sorted(bonds, key=lambda b: min(y.serialNumber for y in b.atoms))

		for b in bonds:
			if b in existing_bondrots_bonds:
				continue
			a = b.atoms[0]
			if  a not in self.nonrotatable and \
				a.idatmType in ('C3', 'N3', 'C2', 'N2') and \
				(a.numBonds > 1 and b.otherAtom(a).numBonds > 1) or \
				a.name == 'DUM' or b.otherAtom(a).name == 'DUM' :
					try:
						br = chimera.BondRot(b)
					except (chimera.error, ValueError), v:
						if "cycle" in str(v):
							continue #discard bonds in cycles!
						else:
							raise
					else:
						br.rotanchor = box.find_nearest(self.donor, b.atoms)
						yield br

	def update_rotatable_bonds(self):
		self.rotatable_bonds = list(self.get_rotatable_bonds())

def enable(**kwargs):
	return Torsion(**kwargs)