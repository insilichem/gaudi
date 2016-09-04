#!/usr/bin/env python
# -*- coding: utf-8 -*-


def test_individual(individual):
    assert individual.genes['Molecule'].compound.mol.numAtoms == 2489

def test_expressed_individual(expressed_individual):
    assert expressed_individual.expressed is True