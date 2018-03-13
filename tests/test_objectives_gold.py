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

import sys
import pytest
from conftest import datapath, expressed
from gaudi.objectives.gold import Gold
from gaudi.genes.molecule import Molecule
import os

if 'CCDC_LICENSE_FILE' not in os.environ:
    os.environ['CCDC_LICENSE_FILE'] = "/QFsoft/applic/CCDC/GOLD5.2/ccdc_licence.dat"
os.environ['PATH'] = "/QFsoft/applic/CCDC/GOLD5.2/bin/:/QFsoft/applic/CCDC/GOLD5.2/GOLD/bin/:" + os.environ.get('PATH', '')


@pytest.mark.skipif('CCDC_LICENSE_FILE' not in os.environ,
                    reason='CCDC GOLD not installed or licensed. Set CCDC_LICENSE_FILE.')
@pytest.mark.parametrize("protein, ligand, affinity", [
    ('5er1_protein.mol2', '5er1_ligand.mol2', 14.58),
])
def test_gold(individual, protein, ligand, affinity):
    individual.genes['Protein'] = Molecule(parent=individual, path=datapath(protein))
    individual.genes['Ligand'] = Molecule(parent=individual, path=datapath(ligand))
    individual.__ready__()
    individual.__expression_hooks__()
    objective = Gold(weight=1.0)
    with expressed(individual):
        assert affinity == objective.evaluate(individual)