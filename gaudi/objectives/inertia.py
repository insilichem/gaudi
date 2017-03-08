#!/usr/bin/env python
# -*- coding: utf-8 -*-

##############
# GaudiMM: Genetic Algorithms with Unrestricted
# Descriptors for Intuitive Molecular Modeling
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
This objective calculates the alignment between the axes of inertia
of the given molecules.
"""

# 3rd party
import numpy as np
# Chimera
from Molecule import atom_positions
# GAUDI
from gaudi.objectives import ObjectiveProvider
from gaudi import parse

def enable(**kwargs):
    kwargs = AxesOfInertia.validate(kwargs)
    return AxesOfInertia(**kwargs)


class AxesOfInertia(ObjectiveProvider):

    """
    Calculates the axes of inertia of given molecules and returns
    their alignment deviation.
    """
    _validate = {
        parse.Required('reference'): parse.Molecule_name,
        parse.Required('targets'): [parse.Molecule_name],
        'threshold': parse.All(parse.Coerce(float), parse.Range(min=0, max=1)),
        'only_primaries': parse.Coerce(bool),
        }
    
    def __init__(self, reference=None, targets=None, only_primaries=False,
                 threshold=0.84, *args, **kwargs):
        ObjectiveProvider.__init__(self, **kwargs)
        self.threshold = threshold
        self._reference = reference
        self._targets = targets

    def reference(self, individual):
        """
        The reference molecule. Usually, the biggest in size
        """
        return individual.find_molecule(self._reference).compound.mol

    def targets(self, individual):
        return [individual.find_molecule(name).compound.mol for name in self._targets]

    def evaluate(self, individual):
        reference = self.reference(individual)
        targets = self.targets(individual)

        all_axes = []
        for target in [reference] + targets:
            axes = calculate_axes_of_inertia(target)
            all_axes.append(axes)

        best_cosines = list(calculate_alignment(all_axes[0], *all_axes[1:]))
        return abs(self.threshold - np.mean(best_cosines))


def calculate_alignment(reference_axis, *probes_axes):
    for primary, secondary, tertiary in probes_axes:
        primary_norm = np.linalg.norm(primary)
        all_cosines = []
        for axis in reference_axis:
            cosine = np.dot(primary, axis) / (primary_norm * np.linalg.norm(axis))
            all_cosines.append(np.abs(cosine))
        yield max(all_cosines)

def calculate_axes_of_inertia(molecule):
    coordinates = atom_positions(molecule.atoms, molecule.openState.xform)
    masses = np.fromiter((a.element.mass for a in molecule.atoms), 
                         dtype='float32', count=molecule.numAtoms)
    inertial_matrix = calculate_inertial_matrix(coordinates, masses).reshape((3,3))   
    eigenvals, eigenvecs = np.linalg.eig(inertial_matrix)
    return eigenvecs[:, eigenvals.argsort()]

def calculate_inertial_matrix(coordinates, masses):
    corrected = coordinates - centroid(coordinates, masses)
    return np.sum(_inertial_matrix_one_atom(c, m) for (c, m) in zip(corrected, masses))

def centroid(coordinates, masses):
    return np.average(coordinates, axis=0, weights=masses)

def _inertial_matrix_one_atom(coordinate, mass):
    x, y, z = coordinate
    M = mass
    return np.array([M*(y*y + z*z), -1*M*y*x,       -1*M*z*x,
                    -1*M*x*y,        M*(x*x + z*z), -1*M*z*y,
                    -1*M*x*z,       -1*M*y*z,        M*(x*x + y*y)])