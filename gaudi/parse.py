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
This module parses and validatesYAML input files into convenient objects
that allow per-attribute access to configuration parameters.
"""

# Python
from __future__ import print_function
import logging
from importlib import import_module
import os
from collections import namedtuple, Mapping
from string import ascii_letters
from random import choice
# External dependencies
import yaml
from munch import Munch, munchify
from voluptuous import *
from voluptuous.humanize import validate_with_humanized_errors
# Own
import gaudi

logger = logging.getLogger(__name__)

#####################################################################
# Some useful validators for other schemas (genes, objectives, etc)
#####################################################################


def AssertList(*validators, **kwargs):
    """
    Make sure the value is contained in a list
    """
    def fn(values):
        if not isinstance(values, (list, tuple)):
            values = [values]
        return [validator(v) for validator in validators for v in values]
    return fn


def Coordinates(v):
    return All([float], Length(min=3, max=3))(v)


def Importable(v):
    try:
        import_module(v)
    except ImportError as e:
        raise Invalid("({}: {}) '{}'".format(e.__class__.__name__, e, v))
    else:
        return v


def Molecule_name(v):
    """
    Ideal implementation:

    def fn(v):
        valid = [i['name'] for i in items if i['module'] == 'gaudi.genes.molecule']
        if v not in valid:
            raise Invalid("{} is not a valid Molecule name".format(v))
        return v
    return fn

    However, I must figure a way to get the gene list beforehand
    """
    return str(v)


MoleculeAtom = namedtuple("MoleculeAtom", ["molecule", "atom"])
MoleculeResidue = namedtuple("MoleculeResidue", ["molecule", "residue"])
namedtuples = {("molecule", "atom"): MoleculeAtom,
               ("molecule", "residue"): MoleculeResidue}
def Named_spec(*names):
    """
    Assert that str is formatted like "Molecule/123", with Molecule being
    a valid name of a Molecule gene and 123 a positive int or *
    """
    def fn(v):
        try:
            name, i = str(v).split('/')
            name.strip()
            if Molecule_name(name):
                if i in ('*', 'last', 'first', 'donor', 'acceptor', 'flagged'):
                    pass
                elif int(i) > 0:
                    i = int(i)
                return namedtuples[(names)](name, i)
            raise ValueError
        except (ValueError, AttributeError):
            raise Invalid("Expected <Molecule name>/<number, *, first or last> but got {}".format(v))
    return fn


def Degrees(v):
    return All(Any(float, int), Range(min=0, max=360))(v)


def ResidueThreeLetterCode(v):
    return All(str, Length(min=3, max=3))(v)


def RelPathToInputFile(inputpath=None):
    if inputpath is None:
        # inputpath = os.environ.get('GAUDI_INPUT_PATH', '')
        inputpath = getattr(gaudi, '__input_path__', '')

    @wraps(RelPathToInputFile)
    def fn(v):
        return os.path.normpath(os.path.join(inputpath, os.path.expanduser(v)))
    return fn

def ExpandUserPathExists(v):
    p = os.path.expanduser(v)
    if os.path.exists(p):
        return p
    raise ValueError("Path {} does not exist".format(p))


def MakeDir(validator):
    def fn(v):
        v = os.path.expanduser(v)
        try:
            os.makedirs(v)
        except OSError:
            if os.path.isfile(v):
                raise Invalid("Could not create directory. '{}' is a file".format(v))
            elif not os.path.isdir(v):
                raise Invalid('Could not create directory' + v)
        return validator(v)
    return fn


class Settings(Munch):

    """
    Parses a YAML input file with PyYAML, validates it with voluptuous and builds a
    attribute-accessible dict with Munch.

    Parameters
    ----------
    path : str
        Path to YAML file
    """

    default_values = {
        'output': {
            'path': '.',
            'name': ''.join(choice(ascii_letters) for _ in range(5)),
            'precision': None,
            'compress': True,
            'history': False,
            'pareto': True,
            'verbose': True,
            'check_every': 10
        },
        'ga': {
            'population': 10,
            'generations': 3,
            'mu': 1,
            'lambda_': 3,
            'mut_eta': 5,
            'mut_pb': 0.50,
            'mut_indpb': 0.75,
            'cx_eta': 5,
            'cx_pb': 0.5,
        },
        'similarity': {
            'type': 'gaudi.similarity.rmsd',
            'args': [['Ligand'], 2.5],
            'kwargs': {}
        },
        'genes': [{}],
        'objectives': [{}]
    }

    schema = {
            Required('output'): {
                'path': MakeDir(RelPathToInputFile()),
                'name': All(basestring, Length(min=1, max=255)),
                'precision': All(Coerce(int), Range(min=0, max=6)),
                'compress': Coerce(bool),
                'history': Coerce(bool),
                'pareto': Coerce(bool),
                'verbose': Coerce(bool),
                'check_every': All(Coerce(int), Range(min=0)),
            },
            'ga': {
                'population': All(Coerce(int), Range(min=2)),
                'generations': All(Coerce(int), Range(min=0)),
                'mu': All(Coerce(float), Range(min=0, max=1)),
                'lambda_': All(Coerce(float), Range(min=0)),
                'mut_eta': All(Coerce(int), Range(min=0)),
                'mut_pb': All(Coerce(float), Range(min=0, max=1)),
                'mut_indpb': All(Coerce(float), Range(min=0, max=1)),
                'cx_eta': All(Coerce(int), Range(min=0)),
                'cx_pb': All(Coerce(float), Range(min=0, max=1)),
            },
            Required('similarity'): {
                'module': basestring,
                'args': list,
                'kwargs': dict
            },
            Required('genes'): All(Length(min=1), [dict]),
            Required('objectives'): All(Length(min=1), [dict])
        }

    def __init__(self, path=None, validation=True):
        data = self.default_values.copy()
        if path is not None:
            gaudi.__input_path__ = os.environ['GAUDI_INPUT_PATH'] = os.path.dirname(path)
            with open(path) as f:
                loaded = yaml.load(f)
            if validation:
                loaded = self.validate(loaded)
            data = deep_update(data, loaded)
        self.update(munchify(data))
        self._path = path

    @property
    def weights(self):
        return [obj.weight for obj in self.objectives]

    @property
    def name_objectives(self):
        return [obj.name for obj in self.objectives]

    def validate(self, data):
        return validate(self.schema, data)

def deep_update(source, overrides):
    """Update a nested dictionary or similar mapping.

    Modify ``source`` in place.
    """
    for key, value in overrides.iteritems():
        if isinstance(value, Mapping) and value:
            returned = deep_update(source.get(key, {}), value)
            source[key] = returned  
        else:
            source[key] = overrides[key]
    return source


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


def validate(schema, data):
    # TODO - print valid keys when not allowed key is used
    try:
        return validate_with_humanized_errors(data, Schema(schema))
    except Exception as e:
        raise e