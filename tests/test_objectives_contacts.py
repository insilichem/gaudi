#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
from conftest import datapath, expressed
from random import random


@pytest.mark.parametrize("path, which, score", [
    ('5dfr_minimized.pdb', 'clashes', 6.7255),
    ('1amb.pdb', 'hydrophobic', -48.9072)
])
def test_interactions(individual, path, which, score):
    from gaudi.objectives.contacts import Contacts
    from gaudi.genes.molecule import Molecule
    individual.genes['Molecule'] = Molecule(parent=individual, path=datapath(path))
    individual.__ready__()
    individual.__expression_hooks__()
    kwargs = dict(probes=['Molecule'], radius=10*random(), which=which, score=score)
    objective = Contacts(**kwargs)
    with expressed(individual):
        assert abs(score - objective.evaluate(individual)) < 0.0001