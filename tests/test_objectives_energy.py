#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
from conftest import datapath, expressed
from gaudi.objectives.energy import Energy
from gaudi.genes.molecule import Molecule


@pytest.mark.parametrize("path, atoms, bonds, residues", [
    ('5dfr_minimized.pdb', 2489, 2523, 159),
    ('1amb.pdb', 351, 357, 22)
])
def test_topology(individual, path, atoms, bonds, residues):
    individual.genes['Molecule'] = Molecule(parent=individual, path=datapath(path))
    individual.__ready__()
    individual.__expression_hooks__()
    objective = Energy()
    with expressed(individual):
        molecules = objective.molecules(individual)
        topology = objective.chimera_molecule_to_openmm_topology(*molecules)
        assert topology.getNumAtoms() == sum(m.numAtoms for m in molecules)
        assert len(list(topology.bonds())) == sum(m.numBonds for m in molecules)
        assert topology.getNumResidues() == sum(m.numResidues for m in molecules)


def energy(individual, path):
    individual.genes['Molecule'] = Molecule(parent=individual, path=datapath(path))
    individual.__ready__()
    individual.__expression_hooks__()
    objective = Energy()
    with expressed(individual):
        return objective.evaluate(individual)


@pytest.mark.parametrize("path, score", [
    ('5dfr_minimized.pdb', 10701.1228976),
    ('1amb.pdb', 4401.90005384),
])
def test_energy(individual, path, score):
    assert (energy(individual, path) - score) < 0.001


@pytest.mark.parametrize("path", [
    '5dfr_minimized.pdb',
    '1amb.pdb',
])
def test_benchmark_energy(benchmark, individual, path):
    benchmark(energy, individual, path)

