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

import os
import pytest
from conftest import datapath, expressed
from gaudi.genes.molecule import Molecule


@pytest.mark.parametrize("path, atoms", [
    ('5dfr_minimized.pdb', 2489),
    ('1amb.pdb', 438),
    ('3pk2_ligand.pdb', 41),
])
def test_molecule(individual, path, atoms):
    absolute_path = os.path.abspath(datapath(path))
    individual.genes['Molecule'] = Molecule(parent=individual, path=absolute_path)
    individual.__ready__()
    individual.__expression_hooks__()
    assert individual.genes['Molecule'].compound.mol.numAtoms == atoms
    assert individual.genes['Molecule'].compound.mol.openedAs[0] == absolute_path
    with expressed(individual):
        assert individual.expressed is True


@pytest.mark.parametrize("protein, ligand", [
    ('5er1_protein.mol2', '5er1_ligand.mol2'),
])
def test_two_molecules(individual, protein, ligand):
    individual.genes['Ligand'] = Molecule(parent=individual, path=datapath(ligand))
    individual.genes['Protein'] = Molecule(parent=individual, path=datapath(protein))
    individual.__ready__()
    individual.__expression_hooks__()
    with expressed(individual):
        assert individual.expressed is True
