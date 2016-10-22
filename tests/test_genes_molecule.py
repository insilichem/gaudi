#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import pytest
from conftest import datapath, expressed
from gaudi.genes.molecule import Molecule


@pytest.mark.parametrize("path, atoms", [
    ('5dfr_minimized.pdb', 2489),
    ('1amb.pdb', 438),
    ('3pk2_ligand.pdb', 41),
])
def test_molecule(individual, path, atoms):
    absolute_path = os.path.abspath(datapath(path))
    individual.genes['Molecule'] = Molecule(parent=individual, path=absolute_path)
    individual.__ready__()
    assert individual.genes['Molecule'].compound.mol.numAtoms == atoms
    assert individual.genes['Molecule'].compound.mol.openedAs[0] == absolute_path
    with expressed(individual):
        assert individual.expressed is True


@pytest.mark.parametrize("protein, ligand", [
    ('5er1_protein.mol2', '5er1_ligand.mol2'),
])
def test_two_molecules(individual, protein, ligand):
    individual.genes['Ligand'] = Molecule(parent=individual, path=datapath('5er1_ligand.mol2'))
    individual.genes['Protein'] = Molecule(parent=individual, path=datapath('5er1_protein.mol2'))
    individual.__ready__()
    with expressed(individual):
        assert individual.expressed is True