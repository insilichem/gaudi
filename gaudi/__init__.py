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
The GAUDI package is comprised of several core modules that establish the base
architecture to build an extensible platform of molecular design.

The main module is :mod:`gaudi.base`, which defines the :class:`gaudi.base.Individual`,
whose instances represent the potential solutions to the proposed problem. Two plugin
packages allow easy customization of how individuals are defined (:package:`genes`) and
how they are evaluated (:package:`objectives`).

:mod:`gaudi.parse` contains parsing utilities to retrieve the configuration files.

:mod:`gaudi.plugin` holds some magic to make the plugin system work.

:mod:`gaudi.box` is a placeholder for several small functions that are used across GAUDI.
"""

from .version import __version__, __version_info__, __author__, __copyright__

# Logging
import logging
try:  # Python 2.7+
    from logging import NullHandler
except ImportError:
    class NullHandler(logging.Handler):

        def emit(self, record):
            pass

logging.getLogger(__name__).addHandler(NullHandler())
