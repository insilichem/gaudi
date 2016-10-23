#!/usr/bin/env python
# -*- coding: utf-8 -*-

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
