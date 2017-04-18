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
This objective is a wrapper around Chimera's `FindHBond`. 
It returns the number of hydrogen bonds that can be formed
between the target molecule and its environment.

.. todo::

    Evaluate the possible HBonds with some kind of function that
    gives a rough idea of the strength (energy) of each of them.

"""

# Python
import logging
# Chimera
import chimera
from FindHBond import findHBonds
from FindHBond.base import filterHBondsBySel
# GAUDI
from gaudi import parse
from gaudi.objectives import ObjectiveProvider
from gaudi.box import draw_interactions

logger = logging.getLogger(__name__)


def enable(**kwargs):
    kwargs = Hbonds.validate(kwargs)
    return Hbonds(**kwargs)


class Hbonds(ObjectiveProvider):

    """
    Hbonds class

    Parameters
    ----------
    probes : list of str
        Names of molecules being object of analysis
    radius : float
        Maximum distance from any point of probe that is searched
        for a possible interaction
    distance_tolerance : float, optional
        Allowed deviation from ideal distance to consider a valid H bond.
    angle_tolerance : float, optional
        Allowed deviation from ideal angle to consider a valid H bond.
    
    Returns
    -------
    int
        Number of detected Hydrogen bonds.
    """

    _validate = {
        parse.Required('probes'): [parse.Molecule_name],
        'radius': parse.All(parse.Coerce(float), parse.Range(min=0)),
        'distance_tolerance': float,
        'angle_tolerance': float
        }
    def __init__(self, probes=None, radius=5.0, distance_tolerance=0.4, angle_tolerance=20.0,
                 *args, **kwargs):
        ObjectiveProvider.__init__(self, **kwargs)
        self._probes = probes
        self.distance_tolerance = distance_tolerance
        self.angle_tolerance = angle_tolerance
        self.radius = radius

    def molecules(self, ind):
        return [m.compound.mol for m in ind._molecules.values()]

    def probes(self, ind):
        return [ind.find_molecule(p).compound.mol for p in self._probes]

    def evaluate(self, ind):
        """
        Find H bonds within self.radius angstroms from self.probes, and return
        only those that interact with probe. Ie, discard those hbonds in that search
        space whose none of their atoms involved are not part of self.probe.
        """
        molecules = self.molecules(ind)
        probe_atoms = [a for m in self.probes(ind) for a in m.atoms]
        test_atoms = self._surrounding_atoms(probe_atoms, molecules)
        hbonds = findHBonds(molecules, cacheDA=self._cache,
                                      donors=test_atoms, acceptors=test_atoms,
                                      distSlop=self.distance_tolerance,
                                      angleSlop=self.angle_tolerance)
        hbonds = filterHBondsBySel(hbonds, probe_atoms, 'any')

        return len(hbonds)

    def display(self, bonds):
        """
        Mock method to show a graphical depiction of the found H Bonds.
        """
        return draw_interactions(bonds, name=self.name, startCol='00FFFF', endCol='00FFFF')

    ###
    def _surrounding_atoms(self, atoms, molecules):
        self.zone.clear()
        self.zone.add(atoms)
        self.zone.merge(chimera.selection.REPLACE,
                        chimera.specifier.zone(self.zone, 'atom', None, self.radius, molecules))
        return self.zone.atoms()
