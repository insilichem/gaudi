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
from gaudi.genes.search import Search, IDENTITY


def search_gene(individual, path, **kwargs):
    individual.genes['Molecule'] = Molecule(parent=individual, path=datapath(path))
    individual.genes['Search'] = search = Search.with_validation(parent=individual,
                                                                 precision=3, **kwargs)
    individual.__ready__()
    individual.__expression_hooks__()
    return search


@pytest.mark.parametrize("path, target, radius, rotate", [
    ('5er1_ligand.mol2', 'Molecule/1', 3.0, False),
    ('5er1_ligand.mol2', 'Molecule/1', 15.0, True),
])
def test_search_translation(individual, path, target, radius, rotate):
    search = search_gene(individual, path, target=target, radius=radius, rotate=rotate)
    target_atom = individual.find_molecule(search.target.molecule).find_atom(search.target.atom)
    original_position = target_atom.xformCoord()
    with expressed(individual):
        expressed_position = target_atom.xformCoord()
        assert original_position.distance(expressed_position) <= radius


# @pytest.mark.parametrize("path,", [
#     ('3pk2_ligand.pdb',),
# ])
# def test_benchmark_search(benchmark, individual, path):
#     @benchmark
#     def run():
#         search_gene(individual, path, radius=5)
#         with expressed(individual):
#             pass
