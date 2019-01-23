#!/usr/bin/env python
# -*- coding: utf-8 -*-

##############
# GaudiMM: Genetic Algorithms with Unrestricted
# Descriptors for Intuitive Molecular Modeling
#
# https://github.com/insilichem/gaudi
#
# Copyright 2019 Jaime Rodriguez-Guerra, Jean-Didier Marechal
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
import numpy as np
from conftest import datapath, expressed
from gaudi.genes.molecule import Molecule
from gaudi.genes.normalmodes import NormalModes
from gaudi.similarity import _rmsd_squared as sqrmsd


def nma_gene(individual, path, **kwargs):
    individual.genes['Molecule'] = Molecule(parent=individual, path=datapath(path))
    individual.genes['NMA'] = nma = NormalModes.with_validation(parent=individual,
                                                                precision=3, **kwargs)
    individual.__ready__()
    individual.__expression_hooks__()
    return nma


@pytest.mark.parametrize("path, target, group_by", [
    ('1amb.pdb', 'Molecule', None),
    ('1amb.pdb', 'Molecule', 'residues'),
    ('1amb.pdb', 'Molecule', 'mass'),
    # ('1amb.pdb', 'Molecule', 'calpha'),
])
def test_normalmodes_prody(individual, path, target, group_by):
    nma = nma_gene(individual, path, target=target, group_by=group_by, n_samples=500)
    observed_rmsds = []
    for sample in nma.NORMAL_MODES_SAMPLES:
        nma.allele = sample
        with expressed(individual):
            observed_rmsd = sqrmsd(nma._original_coords,
                                   nma.molecule.atomCoordinatesArray())
            assert 0 < observed_rmsd
            observed_rmsds.append(observed_rmsd)
    assert np.isclose(np.mean(observed_rmsds), nma.rmsd * nma.rmsd,
                      rtol=0.1, atol=0.001)
