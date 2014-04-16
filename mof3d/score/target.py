#!/usr/bin/python

# MOF3D
# Multi-Objective Force-Field-Free Docking
# A docking module for UCSF Chimera
# Jaime RGP <https://bitbucket.org/jrgp> @ UAB, 2014
import numpy, chimera
def distance(residues, atom_type, target, threshold, wall=True, avg=True):
	if threshold == 'coordinate':
		pass
	elif threshold == 'covalent':
		threshold = chimera.Element.bondLength(chimera.Element(atom_type), target.element)
		wall = True

	distances = {}
	if wall:
		for r in residues:
			ds = [ _distance(a, target) for a in r.atoms if a.element.name == atom_type ]
			if all([d > threshold for d in ds]):
				distances[r] = min(ds, key = lambda x: abs(x-threshold))
			else:
				distances[r] = 1000
	else:
		for r in residues:
			ds = [ _distance(a, target) for a in r.atoms if a.element.name == atom_type ]
			distances[r] = min(ds, key = lambda x: abs(x-threshold))

	if avg:
		return numpy.mean(distances.values())
	return distances

## Internal use
def _distance(atom1, atom2):
	return atom1.xformCoord().distance(atom2.xformCoord())