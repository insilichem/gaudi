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

    Metaclass trickery obtained from `Marty Alchin's blog 
    <http://martyalchin.com/2008/jan/10/simple-plugin-framework/>`_
    Each mount point (ie, ``genes`` and ``objectives``), MUST inherit this one.

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

    for p in plugins:
        if p.name not in container:
            plugin_kwargs = kwargs.copy()
            plugin_kwargs.update(p)
            container[p.name] = sys.modules[p.module].enable(**plugin_kwargs)
            logger.debug("Loaded plugin %s", p.module)
        else:
            logger.error("Already loaded a plugin with same name %s", p.name)
    return container
