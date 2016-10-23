#!/usr/bin/env python
# -*- coding: utf-8 -*-

##############
# GAUDIasm: Genetic Algorithms for Universal
# Design Inference and Atomic Scale Modeling
# Authors:  Jaime Rodriguez-Guerra Pedregal
#            <jaime.rodriguezguerra@uab.cat>
#           Jean-Didier Marechal
#            <jeandidier.marechal@uab.cat>
# Web: https://bitbucket.org/insilichem/gaudi
##############

"""
This module collects more meaningful exceptions than builtins.

"""

class AtomsNotFound(Exception):
    pass

class ResiduesNotFound(Exception):
    pass

class MoleculesNotFound(Exception):
    pass

class TooManyAtoms(Exception):
    pass

class TooManyResidues(Exception):
    pass