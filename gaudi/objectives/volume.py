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
This objective calculates the volume occupied by the
requested Molecule gene instance.

.. note::

    Volume is calculated from the ``surfacePiece`` created by
    a new experimental method found in Chimera's ``Surface.gridsurf``.
    This could be used as an objective for SES, instead of Solvation.

"""
# Python
import logging
# 3rd party
import MeasureVolume
import Surface
# GAUDI
from gaudi import parse
from gaudi.objectives import ObjectiveProvider

logger = logging.getLogger(__name__)


def enable(**kwargs):
    return Volume(**kwargs)


class Volume(ObjectiveProvider):

    """
    Volume class

    Parameters
    ----------
    threshold : float or 'auto'
        Final volume to target. If 'auto', it will calculate the sum of
        VdW volumes of all requested atoms in `probes`. (Unimplemented!)
    target : list of str
        Molecule gene name to calculate volume over

    Returns
    -------
    volume: float
        Calculated volume in A³ (or nm³?)
    """

    validate = parse.Schema({
        'target': [parse.Molecule_name],
        'threshold': [float, 'auto']
        })

    def __init__(self, threshold=None, target=None,
                 *args, **kwargs):
        ObjectiveProvider.__init__(self, **kwargs)
        self.threshold = threshold
        self._target = target

    def target(self, ind):
        return ind.genes[self._target].compound.mol

    def evaluate(self, ind):
        molecule = self.target(ind)
        surface = Surface.gridsurf.ses_surface(molecule.atoms)
        volume, area, holes = MeasureVolume.surface_volume_and_area(surface)
        return abs(volume - self.threshold)
