#!/usr/bin/python

# FETRA (F3T2A)
# Force-field free target-driven docking algorithm
# A docking module for UCSF Chimera
# Jaime RGP <https://bitbucket.org/jrgp> @ UAB, 2014

def distance(atom1, atom2, threshold, wall=True):
	return abs(_distance(atom1, atom2) - threshold)


## Internal use
def _distance(atom1, atom2):
	return atom1.xformCoord().distance(atom2.xformCoord())