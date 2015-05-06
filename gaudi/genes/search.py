#!/usr/bin/python

##############
# GAUDIasm: Genetic Algorithms for Universal
# Design Inference and Atomic Scale Modeling
# Author: Jaime Rodriguez-Guerra Pedregal
# Email: jaime.rogue@gmail.com
# Web: https://bitbucket.org/jrgp/gaudi
##############

from gaudi.genes import GeneProvider
import gaudi.parse
import chimera, random
from chimera import Xform as X
import Matrix as M
from FitMap.search import random_rotation, random_translation_step

ZERO = chimera.Point(0.0, 0.0, 0.0)

def enable(**kwargs):
	return Search(**kwargs)

class Search(GeneProvider):
	
	def __init__(self, target=None, origin=None, radius=None, rotate=True,
				**kwargs):
		GeneProvider.__init__(self, **kwargs)
		self.radius = radius
		self.rotate = rotate
		self.target = target
		mol, serial = gaudi.parse.parse_rawstring(origin)
		try:
			if isinstance(serial, int):
				atom = next(a for a in self.parent.genes[mol].compound.mol.atoms 
							if serial == a.serialNumber)
			else:
				atom = next(a for a in self.parent.genes[mol].compound.mol.atoms 
							if serial == a.name)
		except StopIteration: #atom not found
			raise
		else:
			self.origin = atom.coord().data()

		self.allele = self.random_transform()

	def express(self):
		self.parent.genes[self.target].compound.mol.openState.xform = \
			M.chimera_xform(M.multiply_matrices(*self.allele))
	def unexpress(self):
		self.parent.genes[self.target].compound.mol.openState.xform = \
			X()

	def mate(self, mate):
		xf1 = M.chimera_xform(M.multiply_matrices(*self.allele))
		xf2 = M.chimera_xform(M.multiply_matrices(*mate.allele))
		interp = M.xform_matrix(M.interpolate_xforms(xf1, chimera.Point(0,0,0), 
													xf2, 0.5))
		interp_rot = [ x[:3]+(0,) for x in interp ]
		interp_tl = [ y[:3]+x[-1:] for x, y in zip(interp, M.identity_matrix())]
		self.allele, mate.allele =  (self.allele[0], interp_rot, self.allele[-1]), \
								(interp_tl, mate.allele[1], mate.allele[-1])

	def mutate(self, indpb):
		if random.random() < self.indpb:
			# Careful! Mutation generates a whole NEW position (similar to eta ~= 0)
			# TODO: We could use a eta param in mutation by interpolating original and 
			# a new random xform with a given `frac` parameter
			self.allele = self.random_transform()

	def write(self, path, name):
		pass

	#####
	def random_transform(self):
		ctf = M.translation_matrix([-x for x in self.origin]).tolist()
		rot = random_rotation() if self.rotate else X()
		shift = random_translation_step(self.origin, self.radius).tolist()
		return shift, rot, ctf



### Some useful functions
def translate(molecule, anchor, target):
	if isinstance(anchor, chimera.Atom):
		anchor = anchor.coord()
	if isinstance(target, chimera.Atom):
		target = target.coord()
	
	t = X.translation(target - anchor)
	for a in molecule.atoms:
		a.setCoord(t.apply(a.coord()))

def rotate(molecule, at, alpha):
	if len(at) == 3:
		try:
			a1, a2, a3 = [ a.coord() for a in at ]
		except AttributeError:
			a1, a2, a3 = at
		axis_a = a1 - a2
		axis_b = a3 - a2
		delta = chimera.angle(a1, a2, a3) - alpha
		axis = chimera.cross(axis_a, axis_b)
		if axis.data() == (0.0, 0.0, 0.0):
			axis = chimera.cross(axis_a, axis_b + chimera.Vector(1,0,0))
			print "Warning, had to choose arbitrary normal vector"
		pivot = a2
	elif len(at) == 4:
		try:
			a1, a2, a3, a4 = [ a.coord() for a in at ]
		except AttributeError:
			a1, a2, a3, a4 = at
		axis = a3 - a2
		delta = chimera.dihedral(a1, a2, a3, a4) - alpha
		pivot = a3
	else:
		raise ValueError("Atom list must contain 3 (angle) or 4 (dihedral) atoms only")

	r = X.translation(pivot - ZERO) # move to origin
	r.multiply(X.rotation(axis, - delta)) # rotate
	r.multiply(X.translation(ZERO - pivot)) # return to orig pos
	for a in molecule.atoms:
		a.setCoord(r.apply(a.coord()))

def rand_xform(center, r, rotate=True):
	ctf = M.translation_matrix([-x for x in center]).tolist()
	rot = random_rotation() #if rotate else chimera.Xform().
	shift = random_translation_step(center, r).tolist()
	return shift, rot, ctf