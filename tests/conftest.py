#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pychimera import patch_environ, enable_chimera
patch_environ()
enable_chimera()


import os
import pytest
from gaudi.base import expressed, MolecularIndividual, Environment
from gaudi.box import suppress_ksdssp
import chimera

TESTPATH = os.path.dirname(os.path.abspath(__file__))
chimera.triggers.addHandler("Model", suppress_ksdssp, None)


def datapath(path):
    return os.path.join(TESTPATH, 'data', path)


@pytest.fixture
def individual():
    individual = MolecularIndividual(dummy=True)
    yield individual
    individual.clear_cache()
    chimera.closeSession()

individual2 = individual

@pytest.fixture
def environment():
    environment = Environment()
    yield environment
    environment.clear_cache()
