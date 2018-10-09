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
Objectives are modules that reside in the `gaudi.objectives` package,
and have a certain class structure
"""

# Python
import abc
import logging
from uuid import uuid4
# Chimera
import chimera
# GAUDI
from gaudi import plugin, parse
from gaudi.base import Environment

logger = logging.getLogger(__name__)


class ObjectiveProvider(object):

    """
    Base class that every `objectives` plugin MUST inherit.

    Mount point for plugins implementing new objectives to be evaluated by DEAP.
    The objective resides within the Fitness attribute of the individual.
    Do whatever you want, but use an evaluate() function to return the results.
    Apart from that, there's no requirements.

    The base class includes some useful attributes, so don't forget to call
    `ObjectiveProvider.__init__` in your overriden `__init__`. For example,
    `self.zone` is a `Chimera.selection.ItemizedSelection` object which is shared among
    all objectives. Use that to get atoms in the surrounding of the target gene, and
    remember to `self.zone.clear()` it before use.

    ---
    From (M.A. itself)[http://martyalchin.com/2008/jan/10/simple-plugin-framework/]:
    Now that we have a mount point, we can start stacking plugins onto it.
    As mentioned above, individual plugins will subclass the mount point.
    Because that also means inheriting the metaclass, the act of subclassing
    alone will suffice as plugin registration. Of course, the goal is to have
    plugins actually do something, so there would be more to it than just
    defining a base class, but the point is that the entire contents of the
    class declaration can be specific to the plugin being written. The plugin
    framework itself has absolutely no expectation for how you build the class,
    allowing maximum flexibility. Duck typing at its finest.
    """

    __metaclass__ = plugin.PluginMount
    _cache = {}
    _validate = {}
    _schema = {parse.Required('environment'): Environment,
               'module': parse.Importable,
               'name': str,
               'weight': parse.Coerce(float),
               'zone': chimera.selection.ItemizedSelection,
               'precision': parse.All(parse.Coerce(int), parse.Range(min=0, max=9))}

    def __init__(self, environment=None, name=None, weight=None, zone=None,
                 precision=3, **kwargs):
        self.environment = environment
        self.name = name if name is not None else str(uuid4())
        self.weight = weight
        self.zone = zone if zone is not None else chimera.selection.ItemizedSelection()
        self.precision = precision

    def __ready__(self):
        pass

    @abc.abstractmethod
    def evaluate(self, individual):
        """
        Return the score of the individual under the current conditions.
        """

    @classmethod
    def clear_cache(cls):
        cls._cache.clear()

    @classmethod
    def validate(cls, data, schema=None):
        schema = cls._schema.copy() if schema is None else schema
        schema.update(cls._validate)
        return parse.validate(schema, data)

    @classmethod
    def with_validation(cls, **kwargs):
        cls.__init__(**cls.validate(kwargs))
