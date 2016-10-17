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
The GAUDI package is comprised of several core modules that establish the base
architecture to build an extensible platform of molecular design.

The main module is :mod:`gaudi.base`, which defines the :class:`gaudi.base.Individual`,
whose instances represent the potential solutions to the proposed problem. Two plugin
packages allow easy customization of how individuals are defined (:mod:`gaudi.genes`) and
how they are evaluated (:mod:`gaudi.objectives`).

:mod:`gaudi.parse` contains parsing utilities to retrieve the configuration files.

:mod:`gaudi.plugin` holds some magic to make the plugin system work.

:mod:`gaudi.box` is a placeholder for several small functions that are used across GAUDI.
"""

# Logging
import logging
try:  # Python 2.7+
    from logging import NullHandler
except ImportError:
    class NullHandler(logging.Handler):

        def emit(self, record):
            pass

logging.getLogger(__name__).addHandler(NullHandler())

__author__ = 'Jaime Rodriguez-Guerra, and Jean-Didier Marechal'
__copyright__ = '2016, InsiliChem'
__url__ = 'https://bitbucket.org/insilichem/gaudi'
__description__ = 'GAUDI: Genetic Algorithms for Universal Design Inference'

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions
