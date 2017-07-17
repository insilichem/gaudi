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
from gaudi.genes.molecule import Molecule
from gaudi.genes.rotamers import Rotamers

def rotamers(individual, path, position, seed):
    individual.genes['Molecule'] = Molecule(parent=individual, path=datapath(path))
    individual.genes['Rotamers'] = rotamers = Rotamers(parent=individual,
                                                       residues=[('Molecule', position)])
    individual.__ready__()
    individual.__expression_hooks__()
    rotamers.allele = [seed]
    return rotamers


@pytest.mark.parametrize("path, position, seed, restype, original_chis, new_chis", [
    ('4c3w_protein.mol2', 5, 0, 'ARG', [179.734, 178.061, 60.608, 90.076], 
                                       [179.734, 178.061, 60.608, 90.076]),
    ('4c3w_protein.mol2', 5, 0.015, 'ARG', [179.734, 178.061, 60.608, 90.076], 
                                           [-174.7, -174.4, 69.2, 79.1]),
])
def test_rotamers(individual, path, position, seed, restype, original_chis, new_chis):
    rotamer = rotamers(individual, path, position, seed)
    residue = rotamer.residues[('Molecule', position)]
    # Cache initial coordinates for the residue and ALL the alpha carbons in the protein
    alpha_carbons_coords = [a.xformCoord() for a in residue.molecule.atoms if a.name == 'CA']
    residue_coords = [a.xformCoord() for a in residue.atoms]
    
    with expressed(individual):
        assert residue.id.position == position
        assert residue.type == restype
        assert [a.xformCoord() for a in residue.molecule.atoms if a.name == 'CA'] == alpha_carbons_coords
        assert [a.xformCoord() for a in residue.atoms] != residue_coords
        for real, torsion in zip(original_chis, residue._rotamer_torsions):
            assert abs(real - torsion.chi) < 0.001
        for real, computed_modified in zip(new_chis, rotamer.all_chis(residue)):
            assert abs(real - computed_modified) < 0.001

    # Check everything is back in place!
    for original, reverted in zip([a.xformCoord() for a in residue.atoms], residue_coords):
        for a, b in zip(original, reverted):
            assert abs(a - b) < 0.001

    for original, reverted in zip([a.xformCoord() for a in residue.molecule.atoms if a.name == 'CA'], alpha_carbons_coords):
        for a, b in zip(original, reverted):
            assert abs(a - b) < 0.001

    for original, reverted in zip(original_chis, rotamer.all_chis(residue)):
        assert abs(original - reverted) < 0.001


@pytest.mark.parametrize("path, position, seed", [
    ('4c3w_protein.mol2', 5, 0),
])
def test_benchmark_rotamers(benchmark, individual, path, position, seed):
    @benchmark
    def run():
        rotamers(individual, path, position, seed)
        with expressed(individual):
            pass
