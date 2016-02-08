#!/usr/bin/env python
# -*- coding: utf-8 -*-

##############
# GAUDIasm: Genetic Algorithms for Universal
# Design Inference and Atomic Scale Modeling
# Authors:  Jaime Rodriguez-Guerra Pedregal
#            <jaime.rodriguezguerra@uab.cat>
#           Jean-Didier Marechal
#            <jeandidier.marechal@uab.cat>
# Web: https://bitbucket.org/jrgp/gaudi
##############

"""
This module parses YAML input files into convenient objects that allow
per-attribute access to configuration parameters.

.. todo ::

    Use AttrDict or Bunch instead, and deprecate this shitty code :)

"""

# Python
import logging
# External dependencies
import yaml
from munch import Munch
from voluptuous import (Schema, Required, All, Length, Range, Coerce, Extra, REMOVE_EXTRA)

logger = logging.getLogger(__name__)


class Settings(Munch):

    """
    Parses a YAML input file with PyYAML, validates it with voluptuous and builds a
    attribute-accessible dict with Munch.

    Parameters
    ----------
    path : str
        Path to YAML file
    """

    validate = Schema({
        Required('output'): {
            'path': str,
            'name': All(str, Length(min=1, max=255)),
            'precision': All(int, Range(min=0, max=6)),
            'compress': Coerce(bool),
            'history': Coerce(bool),
            'pareto': Coerce(bool),
        },
        'ga': {
            'population': All(Coerce(int), Range(min=2)),
            'generations': All(Coerce(int), Range(min=0)),
            'mu': All(Coerce(float), Range(min=0, max=1)),
            'lambda_': All(Coerce(float), Range(min=0, max=1)),
            'mut_eta': All(Coerce(int), Range(min=0)),
            'mut_pb': All(Coerce(float), Range(min=0, max=1)),
            'mut_indpb': All(Coerce(float), Range(min=0, max=1)),
            'cx_eta': All(Coerce(int), Range(min=0)),
            'cx_pb': All(Coerce(float), Range(min=0, max=1)),
        },
        'similarity': {
            'type': str,
            'args': [],
            'kwargs': {}
        },
        Required('genes'): All([{Extra: object}], Length(min=1)),
        Required('objectives'): All([{Extra: object}], Length(min=1))
    }, extra=REMOVE_EXTRA)

    default_values = {
        'output': {
            'path': '.',
            'name': '',
            'precision': None,
            'compress': True,
            'history': False,
            'pareto': True,
        },
        'ga': {
            'population': 10,
            'generations': 3,
            'mu': 0.75,
            'lambda_': 0.75,
            'mut_eta': 5,
            'mut_pb': 0.10,
            'mut_indpb': 0.05,
            'cx_eta': 5,
            'cx_pb': 0.75,
        },
        'similarity': {
            'type': 'gaudi.similarity.rmsd',
            'args': [None, 2.5],
            'kwargs': {}
        }
    }

    def __init__(self, path=None):
        self.update(self.default_values)
        if path is not None:
            self._path = path
            with open(path) as f:
                raw_dict = yaml.load(f)
        validated = self.validate(raw_dict)
        self.update(Munch(validated))

    def _weights(self):
        return [obj.weight for obj in self.objectives]

    def _objectives(self):
        return [obj.name for obj in self.objectives]


def parse_rawstring(s):
    """
    It parses reference strings contained in some fields.

    These strings contain references to :class:`genes.molecule` instances, and one of its atoms
    """
    molecule, res_or_atom = s.split('/')
    molecule = molecule.strip()
    try:
        res_or_atom = int(res_or_atom)
    except ValueError:
        pass  # is str
    return molecule, res_or_atom
