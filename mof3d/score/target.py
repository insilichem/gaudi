#!/usr/bin/python

# MOF3D
# Multi-Objective Force-Field-Free Docking
# A docking module for UCSF Chimera
# Jaime RGP <https://bitbucket.org/jrgp> @ UAB, 2014
import numpy, chimera
def distance(atoms, target, threshold, wall=True):
	distances = []
	for a in atoms:
		d = _distance(a, target)
		if threshold == 'covalent':
			threshold = chimera.Element.bondLength(a.element, target.element)
		
		if d < threshold and wall:
			distances.append(1000)
		else:
			distances.append(d-threshold)
	
	return numpy.mean(distances)

## Internal use
def _distance(atom1, atom2):
	return atom1.xformCoord().distance(atom2.xformCoord())