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
This module provides the basic functionality for the plugin system
of `genes` and `objectives`.

"""

# Python
from collections import OrderedDict
from importlib import import_module
import abc
import logging
import sys

logger = logging.getLogger(__name__)


class PluginMount(type):

    """
    Base class for plugin mount points.

    Metaclass trickery obtained from
    [Marty Alchin's blog](http://martyalchin.com/2008/jan/10/simple-plugin-framework/)
    Each mount point (ie,`genes` and `objectives`), MUST inherit this one.

    """

    __metaclass__ = abc.ABCMeta

    def __init__(cls, name, bases, attrs):
        if not hasattr(cls, 'plugins'):
            # This branch only executes when processing the mount point itself.
            # So, since this is a new plugin type, not an implementation, this
            # class shouldn't be registered as a plugin. Instead, it sets up a
            # list where plugins can be registered later.
            cls.plugins = []
        else:
            # This must be a plugin implementation, which should be registered.
            # Simply appending it to the list is all that's needed to keep
            # track of it later.
            cls.plugins.append(cls)


# Helper functions
def import_plugins(*pluginlist):
    """
    Import requested modules, only once, when launch.py is called and the
    configuration is parsed successfully.

    Parameters
    ----------
    pluginlist : list of gaudi.parse.Param
        Usually, the genes or objectives list resulting from
        the configuration parsing.

    """
    plugins = []
    names = set(plugin.module for plugin in pluginlist)
    for name in names:
        try:
            module = import_module(name)
        except ImportError:
            logger.exception('Import error while loading plugin %r', name)
            raise
        else:
            logger.debug('Imported plugin %s', name)
            plugins.append(module)
    return plugins


def load_plugins(plugins, container=None, **kwargs):
    """
    Requests an instance of the class that resides in each plugin. For genes, each
    individual has its own instance, but objectives are treated like a singleton. So,
    they are only instantiated once. That's the reason behind usen a mutable container.

    Parameters
    ----------
    plugins : list of `gaudi.parse.Param`
        Modules to load. Each Param must have a `module` attr with a full import path.
    container : dict or dict-like
        If provided, use this container to retain instances across individuals.
    kwargs : 
        Everything else will be passed to the requested plugin instances.

    """
    if container is None:
        container = OrderedDict()

    for plugin in plugins:
        if plugin.name not in container:
            module = plugin.module
            kwargs.update(plugin.__dict__)
            container[plugin.name] = sys.modules[module].enable(**kwargs)
            logger.debug("Loaded plugin %s", module)
        else:
            logger.error("Already loaded a plugin with same name %s",
                         plugin.name)
    return container
