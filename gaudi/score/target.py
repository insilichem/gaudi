#!/usr/bin/python

# gaudi
# Genetic Algorithm for Unified Docking Inference
# A docking module for UCSF Chimera
# Jaime RGP <https://bitbucket.org/jrgp> @ UAB, 2014
import chimera
def distance(atoms, target, threshold, threshold2=0.0):
	distances = []
	for a in atoms:
		d = _distance(a, target)
		if threshold == 'covalent':
			threshold = chimera.Element.bondLength(a.element, target.element)
		
		if d < threshold2:
			distances.append(1000)
		else:
			distances.append(d-threshold)
	
	return distances

## Internal use
def _distance(atom1, atom2):
	return atom1.xformCoord().distance(atom2.xformCoord())