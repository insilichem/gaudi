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
from random import random
from gaudi.objectives.contacts import Contacts
from gaudi.genes.molecule import Molecule


def interactions(individual, path, which):
    individual.genes['Molecule'] = Molecule(parent=individual, path=datapath(path))
    individual.__ready__()
    individual.__expression_hooks__()
    kwargs = dict(probes=['Molecule'], radius=10*random(), which=which)
    objective = Contacts(**kwargs)
    with expressed(individual):
        return objective.evaluate(individual)

#------------------------------------------------------------------------------
# Assertion tests


@pytest.mark.parametrize("path, which, score", [
    ('5dfr_minimized.pdb', 'clashes', 6.7255),
    ('1amb.pdb', 'hydrophobic', -48.9072)
])
def test_interactions(individual, path, which, score):
    assert abs(score - interactions(individual, path, which)) < 0.0001

#------------------------------------------------------------------------------
# Benchmarking tests


@pytest.mark.parametrize("path, which", [
    ('1amb.pdb', 'clashes'),
    ('1amb.pdb', 'hydrophobic')
])
def test_benchmark_interactions(benchmark, individual, path, which):
    benchmark(interactions, individual, path, which)
