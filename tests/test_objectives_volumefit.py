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
from gaudi.objectives.volumefit import VolumeFit
from gaudi.genes.molecule import Molecule


@pytest.mark.parametrize("molecule, volume, atoms_outside", [
    ('1amb.pdb', '1amb.mrc', 112),
    ('5dfr_minimized.pdb', '1amb.pdb', 1463),
])
def test_volumefit(individual, molecule, volume, atoms_outside):
    individual.genes['Molecule'] = Molecule(parent=individual, path=datapath(molecule))
    individual.__ready__()
    individual.__expression_hooks__()
    objective = VolumeFit(weight=-1.0, probe='Molecule', volume=datapath(volume))
    with expressed(individual):
        assert atoms_outside == objective.evaluate(individual)
