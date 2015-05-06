#!/usr/bin/python

##############
# GAUDIasm: Genetic Algorithms for Universal
# Design Inference and Atomic Scale Modeling
# Author: Jaime Rodriguez-Guerra Pedregal
# Email: jaime.rogue@gmail.com
# Web: https://bitbucket.org/jrgp/gaudi
# PluginMount is based on https://djangosnippets.org/snippets/542/
##############

from collections import OrderedDict
from importlib import import_module
import logging
import sys

class PluginMount(type):
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
def import_plugins(*args, **kwargs):
	plugins = []
	names = set(plugin.type for plugin in args)
	for name in names:
		try:
			module = import_module(name)
		except ImportError:
			logger = logging.getLogger('gaudi')
			logger.warning('Ignoring exception while loading the %r plug-in.', name)
			print '!!Ignoring exception while loading the {} plug-in.'.format(name)
			raise
		else:
			plugins.append(module)
	return plugins

def load_plugins(plugins, container=None, **kwargs):
	if container is None:
		container = OrderedDict()
	
	for plugin in plugins:
		module = plugin.type
		kwargs.update(plugin.__dict__)
		container[plugin.name] = sys.modules[module].enable(**kwargs)
	return container