#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
from conftest import datapath, expressed


@pytest.mark.parametrize("path, energy", [
    ('5dfr_minimized.pdb', 10701.1228976),
    ('1amb.pdb', 4401.90005384),
])
def test_energy(individual, path, energy):
    from gaudi.objectives.energy import Energy
    from gaudi.genes.molecule import Molecule
    individual.genes['Molecule'] = Molecule(path=datapath(path))
    individual.__ready__()
    objective = Energy()
    with expressed(individual):
        result = objective.evaluate(individual)
    assert (result - energy) < 0.001
