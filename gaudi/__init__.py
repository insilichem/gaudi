#!/usr/bin/env python
# -*- coding: utf-8 -*-

##############
# GaudiMM: Genetic Algorithms with Unrestricted
# Descriptors for Intuitive Molecular Modeling
# 
# https://github.com/insilichem/gaudi
#
# Copyright 2017 Jaime Rodriguez-Guerra, Jean-Didier Marechal
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
##############

"""
The GAUDI package is comprised of several core modules that establish the base
architecture to build an extensible platform of molecular design.

The main module is :mod:`gaudi.base`, which defines the :class:`gaudi.base.Individual`,
whose instances represent the potential solutions to the proposed problem. Two plugin
packages allow easy customization of how individuals are defined (:mod:`gaudi.genes`) and
how they are evaluated (:mod:`gaudi.objectives`). Additionally:

- :mod:`gaudi.algorithms` is the place to look for the actual GA implementation
- :mod:`gaudi.box` is a placeholder for several small functions that are used across GaudiMM.
- :mod:`gaudi.parse` contains parsing utilities to retrieve the configuration files.
- :mod:`gaudi.plugin` holds some magic to make the plugin system work.
- :mod:`gaudi.similarity` defines the diversity enhancers.
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
__copyright__ = '2017, InsiliChem'
__url__ = 'https://github.com/insilichem/gaudi'
__title__ = 'GaudiMM'
__description__ = ('GaudiMM: Genetic Algorithms with Unrestricted Descriptors '
                   'for Intuitive Molecular Modeling')

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions
