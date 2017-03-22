#!/usr/bin/env python
# -*- coding: utf-8 -*-

##############
# GaudiMM: Genetic Algorithms with Unrestricted
# Descriptors for Intuitive Molecular Modeling
# 
# http://bitbucket.org/insilichem/gaudi
# 
# Copyright 2017 Jaime Rodriguez-Guerra, Jean-Didier Marechal
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
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
import Surface.gridsurf
import scipy.spatial
import chimera
# GAUDI
from gaudi import parse
from gaudi.objectives import ObjectiveProvider
from gaudi.box import silent_stdout

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
        If True, evaluate cavities volume creating a convex hull and
        calculating the difference between convex hull volume and
        molecule volume

    Returns
    -------
    float
        Calculated volume in A³
    """

    _validate = {
        'target': parse.Molecule_name,
        'threshold': parse.Any(float, 'auto'),
        'cavities': bool
        }

    def __init__(self, threshold=0.0, target=None, cavities=False,
                 *args, **kwargs):
        ObjectiveProvider.__init__(self, **kwargs)
        self.threshold = threshold
        self._target = target
        self.cavities = cavities
        self.evaluate = self.evaluate_convexhull if self.cavities else self.evaluate_volume

    def target(self, ind):
        return ind.find_molecule(self._target).compound.mol

    def evaluate_volume(self, ind):
        molecule = self.target(ind)
        with silent_stdout():
            surface = Surface.gridsurf.ses_surface(molecule.atoms)
        volume, area, holes = MeasureVolume.surface_volume_and_area(surface)
        chimera.openModels.close([surface])
        return abs(volume - self.threshold)

    def evaluate_convexhull(self, ind):
        molecule = self.target(ind)
        with silent_stdout():
            surface = Surface.gridsurf.ses_surface(molecule.atoms)
        volume, area, holes = MeasureVolume.surface_volume_and_area(surface)
        convex_volume = convexhull_volume(surface)
        chimera.openModels.close([surface])
        return convex_volume - abs(volume - self.threshold)


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
        volume : float
            Convex hull volume

    Notes
    -----
        Some systems may produce small volume blobs, resulting in a number
        of different surface pieces. This should be discussed in the future.
        # points = surface.surfacePieces[0].geometry[0]
        # convexhull = scipy.spatial.ConvexHull(points)
    """
    volume = 0
    
    for sp in surface.surfacePieces:
        points = sp.geometry[0]
        convexhull = scipy.spatial.ConvexHull(points)
        volume += convexhull.volume

    return volume