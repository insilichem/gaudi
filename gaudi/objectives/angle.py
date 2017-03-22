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
This objective calculates the angle formed by three
given atoms (or the dihedral, if four atoms are given) and returns
the absolute difference of that angle and the target value.

"""

# Python
from __future__ import print_function
import math
import logging
# Chimera
import chimera
# GAUDI
from gaudi.objectives import ObjectiveProvider
from gaudi import parse

logger = logging.getLogger(__name__)


def enable(**kwargs):
    kwargs = Angle.validate(kwargs)
    return Angle(**kwargs)


class Angle(ObjectiveProvider):

    """
    Angle class

    Parameters
    ----------
    threshold : float
        Optimum angle
    probes : list of str
        Atoms that make the angle, expressed as a series of 
        <molecule_name>/<serial_number> strings

    Returns
    -------
    delta : float
        Deviation from threshold angle, in degrees
    """

    _validate = {
        parse.Required('probes'): parse.AssertList(parse.Named_spec("molecule", "atom")),
        parse.Required('threshold'): parse.Any(parse.Coerce(float), parse.In(['planar']))
        }

    def __init__(self, threshold=None, probes=None, *args, **kwargs):
        ObjectiveProvider.__init__(self, **kwargs)
        self.threshold = threshold
        self._probes = probes

    def probes(self, ind):
        for probe in self._probes:
            mol, serial = probe
            for atom in ind.find_molecule(mol).find_atoms(serial):
                yield atom

    def evaluate(self, ind):
        atoms_coords = [a.xformCoord() for a in self.probes(ind)]
        try:
            angle = chimera.angle(*atoms_coords)
        except TypeError:  # four atoms, means dihedral
            angle = chimera.dihedral(*atoms_coords)
        
        if self.threshold == 'planar':
            return abs(math.sin(math.radians(angle)))
        return abs(self.threshold - angle.real)



# TODO: Probes get lost if rotamers are applied!
