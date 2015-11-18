#!/usr/bin/env python
# -*- coding: utf-8 -*-

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
This objective calculates SASA and/or SESA for the given system (or region).

.. note::

    This objective depends on MoleculeSurface, which in turn depends on MSMS package.
    MSMS is known to fail quite often with large proteins, so expect this to not work most
    of the times. (More details)[http://www.rbvi.ucsf.edu/chimera/docs/UsersGuide/surfprobs.html].

    Hopefully, UCSF Chimera 2.0 will implement a custom alternative to MSMS which won't have these
    problems. At least, that's what they stated \
    (here)[http://www.cgl.ucsf.edu/pipermail/chimera-users/2013-February/008497.html].

"""

# Python
import logging
# Chimera
import chimera
import Measure
import MoleculeSurface
from MoleculeSurface import Surface_Calculation_Error
# GAUDI
from gaudi.objectives import ObjectiveProvider

logger = logging.getLogger(__name__)


def enable(**kwargs):
    return Solvation(**kwargs)


class Solvation(ObjectiveProvider):

    """
    Solvation class

    Parameters
    ----------
    which : {'ses', 'sas'}
        Type of solvation to measure
    target : str
        Name of the molecule gene being analyzed
    """

    def __init__(self, which='ses', target=None,
                 *args, **kwargs):
        ObjectiveProvider.__init__(self, **kwargs)
        self._target = target
        self.which = which

    def target(self, ind):
        return ind.genes[self._target].compound.mol

    def molecules(self, ind):
        return tuple(m.compound.mol for m in ind.genes.values()
                     if m.__class__.__name__ == "Molecule")

    def evaluate(self, ind):
        molecules = self.molecules(ind)
        target = self.target(ind)
        try:
            atoms, ses, sas = self._solvation(self.zone_atoms(target,
                                                              molecules))
        except Surface_Calculation_Error:
            raise Surface_Calculation_Error(
                'Problem with solvation calc. Read this: '
                'http://www.rbvi.ucsf.edu/chimera/docs/UsersGuide/surfprobs.html')
        else:
            if self.which == 'ses':
                surfaces = ses
            elif self.which == 'sas':
                surfaces = sas
            return sum(s for (a, s) in zip(atoms, surfaces) if a in target.atoms)

    def zone_atoms(self, probe, molecules):
        self.zone.clear()
        self.zone.add(probe.atoms)
        self.zone.merge(chimera.selection.REPLACE,
                        chimera.specifier.zone(
                            self.zone, 'atom',
                            None, self.radius,
                            molecules
                        )
                        )
        return self.zone.atoms()
    ###

    @staticmethod
    def _solvation(atoms):
        """
        The actual wrapper around Chimera's own wrapper of MSMS
        """
        xyzr_data = Measure.measure.atom_xyzr(atoms)
        surfaces = MoleculeSurface.xyzr_surface_geometry(xyzr_data)
        # return atoms, ses, sas
        return atoms, surfaces[3][:, 0], surfaces[3][:, 1]
