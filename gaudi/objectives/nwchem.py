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
This objective is a wrapper around NWChem.

It expects an additional input template with the keyword $MOLECULE,
which will be replaced by the currently expressed molecule(s). See ``TEMPLATE``
for an example, which works as a default template if none is provided.

A ~/.nwchemrc file should be present. If you installed NWChem with our conda recipe,
you will find the file in $CONDA_PREFIX/etc/default.nwchemrc. Copy it to your $HOME.
"""

# Python
import os
import subprocess
import logging
import imp
import re
import sys
from tempfile import _get_default_tempdir as default_tempdir, _get_candidate_names as tempnames
from distutils.spawn import find_executable
# GAUDI
from gaudi.objectives import ObjectiveProvider
from gaudi import parse

logger = logging.getLogger(__name__)


TEMPLATE = """
start $TITLE
title "GaudiMM job for $TITLE"
geometry units angstrom
$MOLECULE
end
basis
* library 6-31g
end
task dft energy
 """

def enable(**kwargs):
    kwargs = NWChem.validate(kwargs)
    return NWChem(**kwargs)


class NWChem(ObjectiveProvider):

    """
    NWChem class

    Parameters
    ----------
    targets : list of str
        Molecule name(s) to be processed with NWChem. Small ones!
    template : str, optional
        NWChem input template (or path to a file with such contents) containing
        a $MOLECULE placeholder to be replaced by the currently expressed
        molecule(s) requested in ``targets``, and optionally, a $TITLE
        placeholder to be replaced by the job name. If not provided, it will
        default to the ``TEMPLATE`` example (single-point dft energy).
    parser : str, optional
        Path to a Python script containing a top-level function called
        `parse_output` which will parse the NWChem output and return a
        float. This replaces the default parser, which looks for the last
        'Total <whatever> energy' value.
    processors : int, optional=None
        Number of physical processors to use with openmpi

    Returns
    -------
    float
        Any numeric value as reported by the `parser` routines. By default,
        last 'Total <whatever> energy' value.
    """
    _validate = {
        parse.Required('targets'): [parse.Molecule_name],
        'template': basestring,
        'parser': parse.ExpandUserPathExists,
        'processors': int
        }

    def __init__(self, template=None, targets=('Ligand',), parser=None,
                 executable=None, basis_library=None, processors=None,
                 *args, **kwargs):
        if kwargs.get('precision', 6) < 6:
            kwargs['precision'] = 6
        ObjectiveProvider.__init__(self, **kwargs)
        self.targets = targets
        self.executable = find_executable('nwchem') if executable is None else executable
        self._nprocessors = processors
        self._mpirun = find_executable('mpirun') if processors is not None else None

        if template is None:
            self.template = TEMPLATE
        elif os.path.isfile(template):
            with open(template) as f:
                self.template = f.read()
        else:
            self.template = template
        self.template = self.template.replace('$TITLE', self.environment.cfg.output.name)
        if parser is not None:
            self.parse_output = imp.load_source('_nwchem_parser', parser).parse_output
        if basis_library is not None and os.path.isdir(basis_library):
            if not basis_library.endswith('/'):
                basis_library += '/'
            os.environ['NWCHEM_BASIS_LIBRARY'] = basis_library
        self._oldworkingdir = os.getcwd()
        self._paths = {}
        if os.name == 'posix' and os.path.exists('/dev/shm'):
            self.tmpdir = '/dev/shm'
        else:
            self.tmpdir = default_tempdir()

    def get_molecule_by_name(self, ind, *names):
        """
        Get a molecule gene instance of individual by its name
        """
        for name in names:
            yield ind.find_molecule(name)

    def evaluate(self, ind):
        """
        Run a subprocess calling DSX binary with provided options,
        and parse the results. Clean tmp files at exit.
        """
        self._tmpfile = os.path.join(self.tmpdir, next(tempnames()))
        os.chdir(self.tmpdir)
        molecules = list(self.get_molecule_by_name(ind, *self.targets))
        nwfile = self.prepare_nwfile(*molecules)
        command = []
        if self._mpirun:
            command.extend([self._mpirun, '-n', str(self._nprocessors)])
        command.extend([self.executable, nwfile])
        try:
            p = subprocess.Popen(command, universal_newlines=True,
                                 stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = p.communicate()
            if stderr:
                logger.warning(stderr)
        except subprocess.CalledProcessError:
            logger.warning("Could not run NWChem with command %s", command)
            return -100000 * self.weight
        else:
            return self.parse_output(stdout)
        finally:
            self.clean()
            os.chdir(self._oldworkingdir)

    def prepare_nwfile(self, *molecules):
        xyz = self.get_xyz(*molecules)
        contents = self.template.replace('$MOLECULE', xyz)
        with open(self._tmpfile + '.nw', 'w') as f:
            f.write(contents)
        return self._tmpfile + '.nw'

    def get_xyz(self, *molecules):
        xyz = []
        for m in molecules:
            xyz.extend(self._xyzlines(m))
        return '\n'.join(xyz)

    def _xyzlines(self, molecule):
        lines = []
        for a in molecule.compound.mol.atoms:
            lines.append('{} {} {} {}'.format(a.element.name, *a.xformCoord().data()))
        return lines

    def parse_output(self, stream):
        result = -100000 * self.weight
        for line in stream.splitlines():
            matches = re.search('Total \w+ energy = *([\d\.-]+)', line)
            if matches:
                result = float(matches.group(1))
        return result

    def clean(self):
        os.remove(self._tmpfile + '.nw')
