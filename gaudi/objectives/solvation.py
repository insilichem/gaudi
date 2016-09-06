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
    (here)[http://www.cgl.ucsf.edu/pipermail/chimera-users/2013-February/008497.html]. You can see
    a beta implementation in :mod:`gaudi.objectives.volume`. It's faster, less prone to errors
    but the output is not in agreement with MSMS'.

"""

# Python
import logging
# Chimera
import chimera
import MoleculeSurface
from MoleculeSurface import Surface_Calculation_Error
import MeasureVolume
import Surface.gridsurf

# GAUDI
from gaudi import parse
from gaudi.objectives import ObjectiveProvider
from gaudi.box import silent_stdout

logger = logging.getLogger(__name__)


def enable(**kwargs):
    kwargs = Solvation.validate(kwargs)
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

    validate = parse.Schema({
        parse.Required('target'): parse.Molecule_name,
        'which': parse.In(['msms_ses', 'msms_sas', 'grid']),
        'threshold': parse.All(parse.Coerce(float), parse.Range(min=0))
        }, extra=parse.ALLOW_EXTRA)

    def __init__(self, which='grid', target=None, threshold=0.0,
                 *args, **kwargs):
        ObjectiveProvider.__init__(self, **kwargs)
        self._target = target
        self.which = which
        self.threshold = threshold
        if which in ('ses', 'sas'):
            self.evaluate = self.evaluate_msms
        else:
            self.evaluate = self.evaluate_grid

    def target(self, ind):
        return ind.genes[self._target].compound.mol

    def molecules(self, ind):
        return tuple(m.compound.mol for m in ind._molecules.values())

    def evaluate_msms(self, ind):
        molecules = self.molecules(ind)
        target = self.target(ind)
        try:
            atoms, ses, sas = self._solvation(self.zone_atoms(target, molecules))
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

    def evaluate_grid(self, ind):
        molecule = self.target(ind)
        with silent_stdout():
            surface = Surface.gridsurf.ses_surface(molecule.atoms)
        volume, area, holes = MeasureVolume.surface_volume_and_area(surface)
        chimera.openModels.close([surface])
        return abs(area - self.threshold)

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
        surfaces = MoleculeSurface.msms_geometry(atoms)
        # return atoms, ses, sas
        return atoms, surfaces[3][:, 0], surfaces[3][:, 1]
