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

import pytest
from conftest import datapath, expressed
from gaudi.objectives.coordination import Coordination
from gaudi.genes.molecule import Molecule
import logging
 
logging.basicConfig(level=logging.DEBUG)

@pytest.mark.parametrize("path, metal, ligands, geometry, score", [
    ('4c3w.mol2', 1005, [1013, 1014, 1021, 1023, 1024, 408], 'octahedron', 2.4696),
])
def test_coordination_uncorrected(individual, path, metal, ligands, geometry, score):
    individual.genes['Molecule'] = Molecule(parent=individual, path=datapath(path))
    individual.__ready__()
    individual.__expression_hooks__()
    objective = Coordination(probe=('Molecule', metal), residues=[('Molecule', '*')],
                             atom_elements=('C N O S'.split()), geometry=geometry,
                             weight=-1.0)
    with expressed(individual):
        sphere = objective.coordination_sphere(individual)[:objective.n_vertices]
        assert len(sphere) == len(ligands)
        assert sorted(ligands) == sorted([a.serialNumber for (d, a) in sphere])
        fitness = objective.evaluate(individual)
        assert abs(fitness - score) < 0.0001

@pytest.mark.parametrize("path, metal, ligands, geometry, score", [
    ('4c3w.mol2', 1005, [1013, 1014, 1021, 1023, 1024, 408], 'octahedron', 3.2405),
])
def test_coordination_corrected(individual, path, metal, ligands, geometry, score):
    individual.genes['Molecule'] = Molecule(parent=individual, path=datapath(path))
    individual.__ready__()
    individual.__expression_hooks__()
    objective = Coordination(probe=('Molecule', metal), residues=[('Molecule', '*')],
                             atom_elements=('C N O S'.split()), geometry=geometry,
                             center_of_mass_correction=True, distance_correction=True,
                             weight=-1.0)
    with expressed(individual):
        sphere = objective.coordination_sphere(individual)[:objective.n_vertices]
        assert len(sphere) == len(ligands)
        assert sorted(ligands) == sorted([a.serialNumber for (d, a) in sphere])
        fitness = objective.evaluate(individual)
        assert abs(fitness - score) < 0.0001
