#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
from conftest import datapath, expressed
from gaudi.genes.molecule import Molecule
from gaudi.genes.rotamers import Rotamers


def rotamers(individual, path, position, seed):
    individual.genes['Molecule'] = Molecule(parent=individual, path=datapath(path))
    individual.genes['Rotamers'] = rotamers = Rotamers(parent=individual,
                                                       residues=[('Molecule', position)])
    individual.__ready__()
    individual.__expression_hooks__()
    rotamers.allele = [seed]
    return rotamers

@pytest.mark.parametrize("path, position, seed, restype, original_chis, new_chis", [
    ('4c3w_protein.mol2', 5, 0, 'ARG', [179.734, 178.061, 60.608, 90.076], [-178.1, 179.9, -178.9, -171.1]),
])
def test_rotamers(individual, path, position, seed, restype, original_chis, new_chis):
    rotamer = rotamers(individual, path, position, seed)
    with expressed(individual):
        residue = rotamer.residues[('Molecule', position)]
        assert residue.id.position == position
        assert residue.type == restype
        for real, torsion in zip(original_chis, residue._rotamer_torsions):
            assert abs(real - torsion.chi) < 0.001
        for real, computed_modified in zip(new_chis, rotamer.all_chis(residue)):
            assert abs(real - computed_modified) < 0.001

    for original, reverted in zip(original_chis, rotamer.all_chis(residue)):
        assert abs(original - reverted) < 0.001

@pytest.mark.parametrize("path, position, seed", [
    ('4c3w_protein.mol2', 5, 0),
])
def test_benchmark_rotamers(benchmark, individual, path, position, seed):
    @benchmark
    def run():
        rotamers(individual, path, position, seed)
        with expressed(individual):
            pass
