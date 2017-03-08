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

import pytest
from os.path import expanduser as expand
from conftest import datapath, expressed
from gaudi.genes.molecule import Molecule
from gaudi.objectives.dsx import DSX


def dsx(individual, protein, ligand, **kwargs):
    individual.genes['Protein'] = Molecule(parent=individual, path=datapath(protein))
    individual.genes['Ligand'] = Molecule(parent=individual, path=datapath(ligand))
    individual.__ready__()
    individual.__expression_hooks__()
    options = dict(
        binary=expand('~/.local/DSX/dsx090_and_hotspotsx061_linux/linux64/dsx_linux_64.lnx'),
        potentials=expand('~/.local/DSX/dsx090_and_hotspotsx061_linux/pdb_pot_0511'),
        terms=[True, False, False, True, False],
        proteins=['Protein'], ligands=['Ligand'],
        sorting=1, cofactor_mode=0, with_metals=False)
    options.update(kwargs)
    objective = DSX(**options)
    with expressed(individual):
        return objective.evaluate(individual)

#------------------------------------------------------------------------------
# Assertion tests


@pytest.mark.parametrize("protein, ligand, score", [
    ('5er1_protein.mol2', '5er1_ligand.mol2', -169.462),
])
def test_dsx(individual, protein, ligand, score):
    assert score == dsx(individual, protein, ligand)


@pytest.mark.parametrize("protein, ligand, score", [
    ('3pk2_protein.pdb', '3pk2_ligand_with_metal.mol2', -161.680),
])
def test_dsx_with_metals(individual, protein, ligand, score):
    assert score == dsx(individual, protein, ligand, with_metals=True)

#------------------------------------------------------------------------------
# Benchmarking tests


@pytest.mark.parametrize("protein, ligand", [
    ('3pk2_protein.pdb', '3pk2_ligand.mol2'),
])
def test_benchmark_dsx(benchmark, individual, protein, ligand):
    benchmark(dsx, individual, protein, ligand)


@pytest.mark.parametrize("protein, ligand", [
    ('3pk2_protein.pdb', '3pk2_ligand_with_metal.mol2'),
])
def test_benchmark_dsx_with_metals(benchmark, individual, protein, ligand):
    benchmark(dsx, individual, protein, ligand, with_metals=True)
