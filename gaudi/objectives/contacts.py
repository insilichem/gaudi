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
This objective provides a wrapper around Chimera's
`DetectClash` that detects clashes and contacts.

Clashes are understood as steric conflicts that increases the energy
of the system. They are evaluated as the sum of volumetric overlapping
of the Van der Waals' spheres of the implied atoms.

Contacts are considered as stabilizing, and they are evaluated with a
Lennard-Jones 12-6 like function.

"""

# Python
from __future__ import print_function
import logging
# Chimera
import chimera
import DetectClash
# GAUDI
from gaudi import parse
from gaudi.objectives import ObjectiveProvider

logger = logging.getLogger(__name__)


def enable(**kwargs):
    kwargs = Contacts.validate(kwargs)
    return Contacts(**kwargs)


class Contacts(ObjectiveProvider):

    """
    Contacts class

    Parameters
    ----------
    probes : str
        Name of molecule gene that is object of contacts analysis
    radius : float
        Maximum distance from any point of probes that is searched
        for possible interactions
    which : {'hydrophobic', 'clashes'}
        Type of interactions to measure
    clash_threshold : float, optional
        Maximum overlap of van-der-Waals spheres that is considered as
        a contact (attractive). If the overlap is greater, it's 
        considered a clash (repulsive)
    hydrophobic_threshold : float, optional
        Maximum overlap for hydrophobic patches.
    hydrophobic_elements : list of str, optional, defaults to [C, S]
        Which elements are allowed to interact in hydrophobic patches
    cutoff : float, optional
        If the overlap volume is greater than this, a penalty is applied. 
        Useful to filter bad solutions.
    bond_separation : int, optional
        Ignore clashes or contacts between atoms within n bonds.

    Returns
    -------
    float
        Lennard-Jones-like energy when `which`=`hydrophobic`,
        and volumetric overlap of VdW spheres in A³ if `which`=`clashes`.
    """
    _validate = {
        parse.Required('probes'): [parse.Molecule_name],
        'radius': parse.All(parse.Coerce(float), parse.Range(min=0)),
        'which': parse.In(['hydrophobic', 'clashes']),
        'clash_threshold': parse.Coerce(float),
        'hydrophobic_threshold': parse.Coerce(float),
        'cutoff': parse.Coerce(float),
        'hydrophobic_elements': [basestring],
        'bond_separation': parse.All(parse.Coerce(int), parse.Range(min=2))
        }
    
    def __init__(self, probes=None, radius=5.0, which='hydrophobic',
                 clash_threshold=0.6, hydrophobic_threshold=-0.4, cutoff=0.0,
                 hydrophobic_elements=('C', 'S'), bond_separation=4, *args, **kwargs):
        ObjectiveProvider.__init__(self, **kwargs)
        self.which = which
        self.radius = radius
        self.clash_threshold = clash_threshold
        self.hydrophobic_threshold = hydrophobic_threshold
        self.cutoff = cutoff
        self.hydrophobic_elements = set(hydrophobic_elements)
        self.bond_separation = bond_separation
        self._probes = probes
        if which == 'hydrophobic':
            self.evaluate = self.evaluate_hydrophobic
            self.threshold = hydrophobic_threshold
        else:
            self.evaluate = self.evaluate_clashes
            self.threshold = clash_threshold

    def molecules(self, ind):
        return [m.compound.mol for m in ind._molecules.values()]
        
    def probes(self, ind):
        return [ind.find_molecule(p).compound.mol for p in self._probes]

    def evaluate_clashes(self, ind):
        positive, negative = self.find_interactions(ind)
        clashscore = sum(abs(vol_overlap) for (a1, a2, overlap, vol_overlap) in negative)
        if self.cutoff and clashscore > self.cutoff:
            clashscore = -1000 * self.weight
        return clashscore
        
    def evaluate_hydrophobic(self, ind):
        positive, negative = self.find_interactions(ind)
        return sum(lj_energy for (a1, a2, overlap, lj_energy) in positive)

    def find_interactions(self, ind):
        atoms = self._surrounding_atoms(ind)
        options = dict(test=atoms, intraRes=True, interSubmodel=True,
                       clashThreshold=self.threshold, assumedMaxVdw=2.1,
                       hbondAllowance=0.2, bondSeparation=self.bond_separation)                                          
        clashes = DetectClash.detectClash(atoms, **options)
        return self._analyze_interactions(clashes)

    def _analyze_interactions(self, clashes):
        """
        Interpret contacts provided by DetectClash.

        Parameters
        ----------
        clashes : dict of dict
            Output of DetectClash. It's a dict of atoms, whose values are dicts.
            These subdictionaries contain all the contacting atoms as keys, and
            the respective overlaping length as values.

        Returns
        -------
        positive : list of list
            Each sublist depict an interaction, with four items: the two involved
            atoms, their distance, and their Lennard-Jones score.
        negative : list of list
            Each sublist depict an interaction, with four items: the two involved
            atoms, their distance, and their volumetric overlap.

        .. note ::
            First, collect atoms that can be involved in hydrophobic interactions.
            Namely, C and S.

            Then, iterate the contacting atoms, getting the distances. For each
            interaction, analyze the distance and, based on the threshold, determine
            if it's attractive or repulsive.

            Attractive interactions are weighted with a Lennard-Jones like function
            (``_lennard_jones``), while repulsive attractions are measured with
            the volumetric overlap of the involved atoms' Van der Waals spheres.

        """
        positive, negative = [], []
        for a1, clash in clashes.items():
            for a2, overlap in clash.items():
                # overlap < clash threshold : can be a hydrophobic interaction
                if overlap <= self.clash_threshold:
                    if (a1.element.name in self.hydrophobic_elements
                        and a2.element.name in self.hydrophobic_elements):
                        lj_energy = self._lennard_jones(a1, a2, overlap)
                        positive.append([a1, a2, overlap, lj_energy])
                # overlap > clash threshold : clash!
                else:
                    volumetric_overlap = self._vdw_vol_overlap(a1, a2, overlap)
                    negative.append([a1, a2, overlap, volumetric_overlap])           
        return positive, negative

    def _surrounding_atoms(self, ind):
        """
        Get atoms in the search zone, based on the molecule and the radius
        """
        self.zone.clear()
        self.zone.add([a for m in self.probes(ind) for a in m.atoms])
        self.zone.merge(chimera.selection.REPLACE,
                        chimera.specifier.zone(self.zone, 'atom', None,
                                               self.radius, self.molecules(ind)))
        return self.zone.atoms()

    @staticmethod
    def _lennard_jones(a1, a2, overlap=None):
        """
        VERY rough approximation of a Lennard-Jones score (12-6).

        Parameters
        ----------
        a1, a2 : chimera.Atom
        overlap : float
            Overlapping radii of involved atoms, as provided
            by DetectClash.

        Notes
        -----
        The usual implementation of a LJ potential is:

            LJ = 4*epsilon*(0.25*((r0/r)**12) - 0.5*((r0/r)**6))

        Two approximations are done:
            - The atoms involves are considered equal, hence the
              distance at which the energy is minimum (r0) is just
              the sum of their radii.
            - Epsilon is always 1.  
        """
        r0 = a1.radius + a2.radius
        if overlap is None:
            distance = a1.xformCoord().distance(a2.xformCoord())
        else:
            distance = r0 - overlap
        x = (r0 / distance)**6
        return (x*x - 2*x)

    @staticmethod
    def _vdw_vol_overlap(a1, a2, overlap=None):
        """
        Volumetric overlap of Van der Waals spheres of atoms.

        Parameters
        ----------
        a1, a2 : chimera.Atom
        overlap : float
            Overlapping sphere segment of involved atoms
        .. note ::
            Adapted from Eran Eyal, Comput Chem 25: 712-724, 2004
        """
        PI = 3.14159265359
        if overlap is None:
            d = a1.xformCoord().distance(a2.xformCoord())
        else:
            d = a1.radius + a2.radius - overlap
        if d == 0:
            return 1000
        h_a, h_b = 0, 0
        if d < (a1.radius + a2.radius):
            h_a = (a2.radius ** 2 - (d - a1.radius) ** 2) / (2 * d)
            h_b = (a1.radius ** 2 - (d - a2.radius) ** 2) / (2 * d)

        return (PI / 3) * ((h_a ** 2) * (3 * a1.radius - h_a) + 
                           (h_b ** 2) * (3 * a2.radius - h_b))
