#!/usr/bin/env python
# -*- coding: utf-8 -*-

# from test_fixtures import individual, expressed_individual        

def test_energy(expressed_individual):
    from gaudi.objectives.energy import Energy
    objective = Energy()
    result = objective.evaluate(expressed_individual)
    assert (result - 10701.122823065081) < 0.01
