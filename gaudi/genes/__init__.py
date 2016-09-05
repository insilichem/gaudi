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
Genes are modules that reside in the `gaudi.genes` package,
and have a certain class structure.
"""

# Python
import abc
import logging
import os
import pprint
# GAUDI
from gaudi import plugin, parse

logger = logging.getLogger(__name__)
pp = pprint.PrettyPrinter(4)


class GeneProvider(object):

    """
    Base class that every `genes` plugin MUST inherit.

    The methods listed here are compulsory for all subclasses, since the
    individual will be using them anyway. If it's not relevant in your plugin,
    just define them with a single `pass` statement. Also, don't forget to call
    `GeneProvider.__init__` in your overriden `__init__` function; it registers
    compulsory

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

    # This sole line is the magic behind the plugin system!
    __metaclass__ = plugin.PluginMount
    _cache = {}

    def __init__(self, parent=None, name=None, cache=None,
                 cxeta=5.0, mteta=5.0, indpb=0.75,
                 **kwargs):
        self.parent = parent
        self.name = name
        self.cxeta = cxeta
        self.mteta = mteta
        self.indpb = indpb

    def __ready__(self):
        pass

    @abc.abstractmethod
    def express(self):
        """
        Compile the gene to an evaluable object.
        """

    @abc.abstractmethod
    def unexpress(self):
        """
        Revert expression.
        """

    @abc.abstractmethod
    def mutate(self):
        """
        Perform a mutation on the gene.
        """

    @abc.abstractmethod
    def mate(self, gene):
        """
        Perform a crossover with another gene of the same kind.
        """

    def write(self, path, name, *args, **kwargs):
        """
        Write results of expression to a file representation.
        """
        fullname = os.path.join(path, '{}_{}.txt'.format(name, self.name))
        with open(fullname, 'w') as f:
            f.write(pp.pformat(self.allele))
        return fullname

    @classmethod
    def clear_cache(cls):
        cls._cache.clear()