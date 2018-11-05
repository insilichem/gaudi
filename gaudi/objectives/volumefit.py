#!/usr/bin/env python
# -*- coding: utf-8 -*-

##############
# GaudiMM: Genetic Algorithms with Unrestricted
# Descriptors for Intuitive Molecular Modeling
#
# https://github.com/insilichem/gaudi
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
This objective calculates how a given molecule fits
within a given volume. Useful for XRC, SAXS, and so on.

"""
# Python
import logging
# 3rd party
import chimera
from FitMap.fitmap import points_outside_contour
# GAUDI
from gaudi import parse
from gaudi.objectives import ObjectiveProvider

logger = logging.getLogger(__name__)


def enable(**kwargs):
    kwargs = VolumeFit.validate(kwargs)
    return VolumeFit(**kwargs)


class VolumeFit(ObjectiveProvider):

    """
    VolumeFit class

    Parameters
    ----------
    probe : str
        Name of the molecule that must be fitted within a volume
    volume : str
        Path to a volume file, or a equivalent molecule (ie SAXS probes)
    resolution : float
        Interpolation resolution (in A) to use for map generation
        if a molecule is provided in ``volume``.

    Returns
    -------
    float
        Number of atoms outside volume (to be minimized)
    """

    _validate = {
        parse.required('probe'): parse.Molecule_name,
        parse.required('volume'): parse.RelPathToInputFile(),
        'resolution': parse.All(parse.Coerce(float), parse.Range(min=0))
        }

    def __init__(self, probe=None, volume=None, resolution=10.0, *args, **kwargs):
        ObjectiveProvider.__init__(self, **kwargs)
        self._volume = volume
        self._probe = probe
        self.resolution = resolution
        self._id_xform = chimera.Xform()

    @property
    def volume(self):
        if not self._cache.get(self.name):
            model = chimera.openModels.open(self._volume)[0]
            if isinstance(model, chimera.Molecule):  # we want the molmap!
                molecule = model
                model = molecule_map(model.atoms, self.resolution)
                chimera.openModels.close(molecule)
            self._cache[self.name] = model
        return self._cache[self.name]

    def probe(self, ind):
        return ind.find_molecule(self._probe).compound.mol

    def evaluate(self, ind):
        return float(points_outside_contour(self.probe(ind).xyz, self._id_xform, self.volume)[0])
