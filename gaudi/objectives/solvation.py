#!/usr/bin/python

##############
# GAUDIasm: Genetic Algorithms for Universal
# Design Inference and Atomic Scale Modeling
# Authors:  Jaime Rodriguez-Guerra Pedregal
#            <jaime.rodriguezguerra@uab.cat>
#           Jean-Didier Marechal
#            <jeandidier.marechal@uab.cat>
# Web: https://bitbucket.org/jrgp/gaudi
##############

"""
Calculates SASA and/or SESA for the given system (or region).

.. note::

    This objective depends on MoleculeSurface, which in turn depends on MSMS package.
    MSMS is known to fail quite often with large proteins, so expect this to not work most
    of the times. (More details)[http://www.rbvi.ucsf.edu/chimera/docs/UsersGuide/surfprobs.html].

    Hopefully, UCSF Chimera 2.0 will implement a custom alternative to MSMS which won't have these
    problems. At least, that's what they stated \
    (here)[http://www.cgl.ucsf.edu/pipermail/chimera-users/2013-February/008497.html].
"""

# Chimera
import Measure
import MoleculeSurface
from MoleculeSurface import Surface_Calculation_Error
# GAUDI
from gaudi.objectives import ObjectiveProvider


def enable(**kwargs):
    return Solvation(**kwargs)


class Solvation(ObjectiveProvider):

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
            return sum(s for (a, s) in zip(atoms, surfaces) if a in self.target.atoms)

    ###
    @staticmethod
    def _solvation(atoms):
        xyzr_data = Measure.measure.atom_xyzr(atoms)
        surfaces = MoleculeSurface.xyzr_surface_geometry(xyzr_data)
        # return atoms, ses, sas
        return atoms, surfaces[3][:, 0], surfaces[3][:, 1]
