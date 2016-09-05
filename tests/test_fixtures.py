#!/usr/bin/env python
# -*- coding: utf-8 -*-

individual_global = None
def test_individual(individual):
    global individual_global
    individual_global = individual
    assert True

def test_individuals(individual):
    global individual_global
    assert id(individual) != id(individual_global)
