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
:mod:`gaudi.cli.gaudi_view` is a wrapper around several GUI programs
that can help visualize GaudiMM results.

As of now, it implements:

- GaudiView: An extension for UCSF Chimera.

"""

from __future__ import print_function
import os
import sys
from subprocess import call
from pychimera.pychimera import guess_chimera_path 


def launch(filename, viewer=None):
    if viewer in (None, 'gaudiview'):
        visualize_with_gaudiview(filename)
    else:
        sys.exit("Viewer {} not supported".format(viewer))


def visualize_with_gaudiview(filename):
    chimera_path = 'chimera'
    chimera_paths = guess_chimera_path(common_locations=True)
    for path in chimera_paths:
        if 'headless' not in path.lower():
            chimera_path = os.path.join(path, 'bin', 'chimera')
    call([chimera_path, filename])
