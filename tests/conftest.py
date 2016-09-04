#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import pytest

TESTPATH = os.path.dirname(os.path.abspath(__file__))


@pytest.fixture
def individual(request):
    from gaudi.base import Individual
    from gaudi.genes.molecule import Molecule
    individual = Individual(dummy=True)
    path = os.path.join(TESTPATH, 'data', '5dfr_minimized.pdb')
    individual.genes['Molecule'] = Molecule(path=path)
    return individual

@pytest.fixture
def expressed_individual(individual):
    individual.express()
    yield individual
    individual.unexpress()

