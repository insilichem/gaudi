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

logger = logging.getLogger(__name__)


class Settings(object):

    """ 
    Simple parser for YAML settings file.

    It should be made of dictionaries and lists of dictionaries. The keys
    will be added as attributes of the returned object.

    Parameters
    ----------
    path : str
        Path to the yaml file.

    """

    def __init__(self, path, asDict=False):
        self._path = path
        self.data = {}
        self._parse()
        self.weights = self._weights()
        self.objectivesnames = self._objectives()

    def _weights(self):
        return [obj.weight for obj in self.objectives]

    def _objectives(self):
        return [obj.name for obj in self.objectives]

    def _parse(self):
        with open(self._path, 'r') as f:
            self.data = yaml.load(f)
        # make dict available as attrs
        for k, v in self.data.items():
            if isinstance(v, list):  # objectives is a list!
                self.__dict__[k] = [Param(d) for d in v]
            else:
                self.__dict__[k] = Param(v)


class Param(object):

    """
    Blank object for storing the attributes used through the parsing

    .. todo::

        Maybe a ``namedtuple`` is better suit for this task?
    """

    def __init__(self, *d):
        for d_ in d:
            self.__dict__.update(d_)


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


#####


def _test_rebuild(cfg):
    for s, c in cfg.items():
        if isinstance(c, dict):
            print '\n[' + s + ']'
            for k, v in c.items():
                print k, "=", v
        elif isinstance(c, list):
            for i, l in enumerate(c):
                print '\n[' + str(s) + ' ' + str(i) + ']'
                for k, v in l.items():
                    print k, "=", v

if __name__ == '__main__':
    import sys
    cfg = Settings(sys.argv[1])
    print[o.module for o in cfg.objective]
