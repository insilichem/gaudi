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
import scipy.spatial
# GAUDI
from gaudi import parse
from gaudi.objectives import ObjectiveProvider

logger = logging.getLogger(__name__)


def enable(**kwargs):
    kwargs = Volume.validate(kwargs)
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
    cavities : boolean, optional, default=False
        if True, evaluate cavities volume creating a convex hull and
        calculating the difference between convex hull volume and
        molecule volume

    Returns
    -------
    volume: float
        Calculated volume in A³ (or nm³?)
    """

    validate = parse.Schema({
        'target': [parse.Molecule_name],
        'threshold': parse.Any(float, 'auto'),
        'cavities': parse.All(bool)
        }, extra=parse.ALLOW_EXTRA)

    def __init__(self, threshold=None, target=None, cavities=False,
                 *args, **kwargs):
        ObjectiveProvider.__init__(self, **kwargs)
        self.threshold = threshold
        self._target = target
        self.cavities = cavities

    def target(self, ind):
        return ind.genes[self._target].compound.mol

    def evaluate(self, ind):
        molecule = self.target(ind)
        surface = Surface.gridsurf.ses_surface(molecule.atoms)
        volume, area, holes = MeasureVolume.surface_volume_and_area(surface)
        if self.cavities:
            return convexhull_volume(surface) - volume
        else:
            return abs(volume - self.threshold)

###
def convexhull_volume(surface):
    """
    This function gets a surface, creates the convex hull and calculates
    its volume

    Parameters
    ----------
        surface : Surface.gridsurf.ses_surface(molecule.atoms)

    Returns
    -------
        convexhull_volume : float
            convex hull volume
    """
    convexhull_volume = 0
    # points = surface.surfacePieces[0].geometry[0]
    # convexhull = scipy.spatial.ConvexHull(points)
    for sp in surface.surfacePieces:
        points = sp.geometry[0]
        convexhull = scipy.spatial.ConvexHull(points)
        convexhull_volume += convexhull.volume
        # sp.geometry = points, convexhull.simplices
    return convexhull_volume