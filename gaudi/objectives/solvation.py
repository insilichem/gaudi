#!/usr/bin/python

##############
# GAUDIasm: Genetic Algorithms for Universal
# Design Inference and Atomic Scale Modeling
# Author: Jaime Rodriguez-Guerra Pedregal
# Email: jaime.rogue@gmail.com
# Web: https://bitbucket.org/jrgp/gaudi
#############

from gaudi.objectives import ObjectiveProvider
import MoleculeSurface, Measure
from MoleculeSurface import Surface_Calculation_Error

def enable(**kwargs):
	return Solvation(**kwargs)

class Solvation(ObjectiveProvider):
	""" This objective depends on MoleculeSurface, which in turn depends on MSMS package.
	MSMS is known to fail quite often with large proteins, so expect this to not work most
	of the times. More details at http://www.rbvi.ucsf.edu/chimera/docs/UsersGuide/surfprobs.html """

	def __init__(self, which='ses', target=None,
				*args, **kwargs):
		ObjectiveProvider.__init__(self, **kwargs)
		self.which = which
		try:
			self.target = self.parent.genes[target].compound.mol
		except KeyError:
			raise
		except AttributeError:
			raise

	def evaluate(self):
		try:
			atoms, ses, sas = self._solvation(self.env.atoms())
		except Surface_Calculation_Error:
			raise Surface_Calculation_Error("""Problem with solvation calc. 
Read this: http://www.rbvi.ucsf.edu/chimera/docs/UsersGuide/surfprobs.html""")
		else:
			if self.which == 'ses':
				surfaces = ses
			elif self.which == 'sas':
				surfaces = sas
			return sum(s for (a,s) in zip(atoms, surfaces) if a in self.target.atoms)

	###
	@staticmethod
	def _solvation(atoms):
		xyzr_data = Measure.measure.atom_xyzr(atoms)
		surfaces = MoleculeSurface.xyzr_surface_geometry(xyzr_data)
		# return atoms, ses, sas
		return atoms, surfaces[3][:,0], surfaces[3][:,1]