#!/usr/bin/python

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
:mod:`gaudi.parse` parses YAML input files into convenient objects that allow
per-attribute access to configuration parameters.
"""

# Python
import os
import logging
# External dependencies
import yaml


class Settings(object):

    """ Simple parser for YAML settings file.

    It should be made of dictionaries and
    lists of dictionaries. The keys will be added as attributes of the returned object.
    :path:  Path to the yaml file.
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


def enable_logging(output_dir=None):
    logger = logging.getLogger('gaudi')
    logger.setLevel(logging.DEBUG)

    # create CONSOLE handler and set level to error
    handler = logging.StreamHandler()
    handler.setLevel(logging.ERROR)
    formatter = logging.Formatter("%(levelname)s - %(message)s")
    handler.setFormatter(formatter)
    logger.addHandler(handler)

    # create debug file handler and set level to debug
    if output_dir:
        handler = logging.FileHandler(
            os.path.join(output_dir, "debug.log"), "w")
        handler.setLevel(logging.DEBUG)
        formatter = logging.Formatter(
            "%(asctime)s %(levelname)s - %(message)s")
        handler.setFormatter(formatter)
        logger.addHandler(handler)

    return logger

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
    print [o.type for o in cfg.objective]
