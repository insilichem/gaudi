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
Objectives are modules that reside in the `gaudi.objectives` package,
and have a certain class structure
"""

# Python
import abc
import logging
from uuid import uuid4
# GAUDI
from gaudi import plugin

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

    def __init__(self, environment=None, name=None, weight=None, zone=None,
                 **kwargs):
        self.environment = environment
        self.name = name if name is not None else str(uuid4())
        self.weight = weight
        self.zone = zone

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
