#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
from conftest import datapath, expressed
from gaudi.objectives.distance import Distance
from gaudi.genes.molecule import Molecule
from gaudi.genes.torsion import Torsion


@pytest.mark.parametrize("path, angle, bonds, rotatable, distance", [
    ('3pk2_ligand.pdb', 1.0, 49, 8, 3.796920067633768),
])
def test_torsion(individual, path, angle, bonds, rotatable, distance):
    individual.genes['Molecule'] = Molecule(parent=individual, path=datapath(path))
    individual.genes['Torsion'] = torsion = Torsion(parent=individual, target='Molecule')
    individual.__ready__()

    torsion.allele = [angle] * torsion.max_bonds

    with expressed(individual):
        assert individual.genes['Molecule'].compound.mol.numBonds == bonds
        atom1 = individual.genes['Molecule'].compound.mol.atoms[0]
        atom2 = individual.genes['Molecule'].compound.mol.atoms[-1]
        assert len(list(torsion.rotatable_bonds)) == rotatable
        assert Distance._distance(atom1, atom2) == distance


@pytest.mark.parametrize("path, angle, distance", [
    ('1amb.pdb', 90.0, 58.986905376603815),
])
def test_backbone_torsion(individual, path, angle, distance):
    individual.genes['Molecule'] = mol = Molecule(parent=individual, path=datapath(path))
    individual.genes['Torsion'] = torsion = Torsion(parent=individual, target='Molecule')
    individual.__ready__()

    torsion.max_bonds = mol.compound.mol.numBonds
    torsion.allele = [angle] * torsion.max_bonds
    torsion.rotatable_atom_types = ()
    torsion.rotatable_atom_names = ('CA',)

    with expressed(individual):
        assert all('CA' in [a.name for a in br.bond.atoms] for br in torsion.rotatable_bonds)
        atom1 = individual.genes['Molecule'].compound.mol.atoms[0]
        atom2 = individual.genes['Molecule'].compound.mol.atoms[-1]
        assert Distance._distance(atom1, atom2) == distance
