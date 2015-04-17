#!/usr/bin/python

# gaudi
# Genetic Algorithm for Unified Docking Inference
# A docking module for UCSF Chimera
# Jaime RGP <https://bitbucket.org/jrgp> @ UAB, 2014
import chimera
def distance(atoms, target, threshold):
	distances = []
	for a in atoms:
		d = _distance(a, target)
		if threshold == 'covalent':
			threshold = chimera.Element.bondLength(a.element, target.element)
		distances.append(d-threshold)
	return distances

def angle():
	pass

def dihedral():
	pass
	
## Internal use
def _distance(atom1, atom2):
	return atom1.xformCoord().distance(atom2.xformCoord())