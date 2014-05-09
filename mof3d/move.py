#!/usr/bin/python

# MOF3D
# Multi-Objective Force-Field-Free Docking
# A docking module for UCSF Chimera
# Jaime RGP <https://bitbucket.org/jrgp> @ UAB, 2014

import chimera
from chimera import Xform as x
from FitMap.search import random_rotation, random_translation_step
from PDBmatrices.matrices import chimera_xform as X
ZERO = chimera.Point(0.0, 0.0, 0.0)

def translate(molecule, anchor, target):
	if isinstance(anchor, chimera.Atom):
		anchor = anchor.coord()
	if isinstance(target, chimera.Atom):
		target = target.coord()
	
	t = x.translation(target - anchor)
	for a in molecule.atoms:
		a.setCoord(t.apply(a.coord()))

def rotate(molecule, at, alpha):
	if len(at) == 3:
		a1, a2, a3 = [ a.coord() for a in at ]
		axis_a = a1 - a2
		axis_b = a3 - a2
		delta = chimera.angle(a1, a2, a3) - alpha
		axis = chimera.cross(axis_a, axis_b)
		if axis.data() == (0.0, 0.0, 0.0):
			axis = chimera.cross(axis_a, axis_b + chimera.Vector(1,0,0))
			print "Warning, had to choose arbitrary normal vector"
		pivot = a2
	elif len(at) == 4:
		a1, a2, a3, a4 = [ a.coord() for a in at ]
		axis = a3 - a2
		delta = chimera.dihedral(a1, a2, a3, a4) - alpha
		pivot = a3
	else:
		raise ValueError("Atom list must contain 3 (angle) or 4 (dihedral) atoms only")

	r = x.translation(pivot - ZERO) # move to origin
	r.multiply(x.rotation(axis, - delta)) # rotate
	r.multiply(x.translation(ZERO - pivot)) # return to orig pos
	for a in molecule.atoms:
		a.setCoord(r.apply(a.coord()))

def rand_xform(mol, center, r):
	has, bbox = mol.bbox()
	centroid = bbox.center()
	xf = x.translation(centroid - ZERO)
	xf.multiply(X(random_rotation()))
	randompos = random_translation_step(center, r)[:,3].tolist()
	xf.translate(ZERO - chimera.Point(*randompos))
	return [xf.getOpenGLMatrix()[i::4] for i in range(3)]
	

