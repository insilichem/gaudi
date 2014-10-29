#!/usr/bin/python

# gaudi
# Genetic Algorithm for Unified Docking Inference
# A docking module for UCSF Chimera
# Jaime RGP <https://bitbucket.org/jrgp> @ UAB, 2014

import deap.tools

def hetCrossover(ind1, ind2):
	for key in ind1:
		if key == 'molecule': 
			continue
		elif key == 'linker_rots':
			ind1[key][:], ind2[key][:] = deap.tools.cxSimulatedBinaryBounded(
				ind1[key], ind2[key], eta=10., low=0., up=360.)
		elif key == 'mutamers':
			ind1[key], ind2[key] = deap.tools.cxTwoPoint(ind1[key], ind2[key])
		elif key == 'rotamers':
			ind1[key], ind2[key] = deap.tools.cxTwoPoint(ind1[key], ind2[key])
	return ind1, ind2

def hetMutation(ind, indpb,):
	for key in ind:
		if key == 'molecule': 
			continue
		elif key == 'linker_rots':
			ind[key] = deap.tools.mutPolynomialBounded(ind[key], 
				eta=10., low=0., up=360., indpb=indpb)[0]
		elif key == 'mutamers':
			ind[key] = deap.tools.mutUniformInt(ind[key], 
				low=0, up=8, indpb=indpb)[0]
		elif key == 'rotamers':
			ind[key] = deap.tools.mutUniformInt(ind[key], 
				low=0, up=8, indpb=indpb)[0]
	return ind,