#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import pytest
from gaudi.base import expressed
import chimera

TESTPATH = os.path.dirname(os.path.abspath(__file__))

def datapath(path):
    return os.path.join(TESTPATH, 'data', path)

@pytest.fixture
def individual():
    from gaudi.base import Individual
    individual = Individual(dummy=True)
    yield individual
    individual.clear_cache()
    chimera.closeSession()

@pytest.fixture
def environment():
    from gaudi.base import Environment
    environment = Environment()
    yield environment
    environment.clear_cache()
