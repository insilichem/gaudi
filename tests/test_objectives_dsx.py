#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
from os.path import expanduser as expand
from conftest import datapath, expressed
from gaudi.genes.molecule import Molecule
from gaudi.objectives.dsx import DSX


@pytest.mark.parametrize("protein, ligand, energy", [
    ('5er1_protein.mol2', '5er1_ligand.mol2', -169.462),
])
def test_dsx(individual, protein, ligand, energy):
    individual.genes['Protein'] = Molecule(parent=individual, path=datapath(protein))
    individual.genes['Ligand'] = Molecule(parent=individual, path=datapath(ligand))
    individual.__ready__()
    individual.__expression_hooks__()
    kwargs = dict(
        binary=expand('~/.local/DSX/dsx090_and_hotspotsx061_linux/linux64/dsx_linux_64.lnx'),
        potentials=expand('~/.local/DSX/dsx090_and_hotspotsx061_linux/pdb_pot_0511'),
        terms=[True, False, False, True, False],
        proteins=['Protein'], ligands=['Ligand'],
        sorting=1, cofactor_mode=0, with_metals=False)
    objective = DSX(**kwargs)
    with expressed(individual):
        assert energy == objective.evaluate(individual)


@pytest.mark.parametrize("protein, ligand, energy", [
    ('3pk2_protein.pdb', '3pk2_ligand_with_metal.mol2', -161.680),
])
def test_dsx_with_metals(individual, protein, ligand, energy):
    individual.genes['Protein'] = Molecule(parent=individual, path=datapath(protein))
    individual.genes['Ligand'] = Molecule(parent=individual, path=datapath(ligand))
    individual.__ready__()
    individual.__expression_hooks__()
    kwargs = dict(
        binary=expand('~/.local/DSX/dsx090_and_hotspotsx061_linux/linux64/dsx_linux_64.lnx'),
        potentials=expand('~/.local/DSX/dsx090_and_hotspotsx061_linux/pdb_pot_0511'),
        terms=[True, False, False, True, False],
        proteins=['Protein'], ligands=['Ligand'],
        sorting=1, cofactor_mode=0, with_metals=True)
    objective = DSX(**kwargs)
    with expressed(individual):
        assert energy == objective.evaluate(individual)