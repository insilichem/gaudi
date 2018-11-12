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
from gaudi.objectives.distance import Distance
from gaudi.genes.molecule import Molecule
from gaudi.genes.torsion import Torsion


def torsions(individual, path, angle, **kwargs):
    individual.genes['Molecule'] = Molecule(parent=individual, path=datapath(path))
    individual.genes['Torsion'] = torsion = Torsion(parent=individual, target='Molecule', **kwargs)
    individual.__ready__()
    individual.__expression_hooks__()
    torsion.allele = [angle] * torsion.max_bonds
    return torsion


@pytest.mark.parametrize("path, angle, bonds, rotatable, distance", [
    ('3pk2_ligand.pdb', 1.0, 49, 8, 3.796920067633768),
])
def test_torsion(individual, path, angle, bonds, rotatable, distance):
    torsion = torsions(individual, path, angle)
    with expressed(individual):
        assert individual.genes['Molecule'].compound.mol.numBonds == bonds
        atom1 = individual.genes['Molecule'].compound.mol.atoms[0]
        point = individual.genes['Molecule'].compound.mol.atoms[-1].xformCoord()
        assert len(list(torsion.rotatable_bonds)) == rotatable == torsion.max_bonds == len(torsion.allele)
        assert Distance._distance(atom1, point) == distance


@pytest.mark.parametrize("path, angle, distance", [
    ('1amb.pdb', 90.0, 58.986905376603815),
])
def test_backbone_torsion(individual, path, angle, distance):
    torsion = torsions(individual, path, angle, rotatable_atom_types=(), rotatable_atom_names=('CA',))
    with expressed(individual):
        assert all('CA' in [a.name for a in br.bond.atoms] for br in torsion.rotatable_bonds)
        atom1 = individual.genes['Molecule'].compound.mol.atoms[0]
        point = individual.genes['Molecule'].compound.mol.atoms[-1].xformCoord()
        assert abs(Distance._distance(atom1, point) - distance) < 0.0001


@pytest.mark.parametrize("path, selection, n_bonds", [
    ('1amb.pdb', '', 82), # n_bonds = 3 * CA - number of gly
    ('1amb.pdb', '@CB,CA', 56), # remove CA-CB bonds
    ('1amb.pdb', '@CB,CA|:0-14', 28), # remove half the protein
])
def test_torsion_backbone_selection(individual, path, selection, n_bonds):
    """ Use `non_rotatable_selection` to discard useless rotations in
        side chains. By selecting '@CB,CA' the bond between them is
        considered and discarded from the set of bonds that contains
        atoms named CA. As a result we only rotate phi and psi angles.
    """
    torsion = torsions(individual, path, 90.0,
                rotatable_atom_types=(),
                rotatable_elements=(),
                rotatable_atom_names=('CA',),
                non_rotatable_selection=selection)
    with expressed(individual):
        assert all('CA' in [a.name for a in br.bond.atoms] for br in torsion.rotatable_bonds)
        assert len(torsion.rotatable_bonds) == n_bonds


@pytest.mark.parametrize("path, angle", [
    ('3pk2_ligand.pdb', 90.0),
])
def test_benchmark_torsion(benchmark, individual, path, angle):
    @benchmark
    def run():
        torsions(individual, path, angle)
        with expressed(individual):
            pass
