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
This objective is a wrapper around the scoring functions provided by
`AutoDock Vina <http://vina.scripps.edu/>`_.
"""

import logging
import os
import sys
from subprocess import call, check_output, CalledProcessError
from tempfile import _get_default_tempdir, _get_candidate_names
from gaudi.objectives import ObjectiveProvider
from gaudi import parse
import MolKit
from AutoDockTools.MoleculePreparation import AD4ReceptorPreparation, AD4LigandPreparation

logger = logging.getLogger(__name__)


def enable(**kwargs):
    kwargs = Vina.validate(kwargs)
    return Vina(**kwargs)


class Vina(ObjectiveProvider):

    """
    Vina class

    Parameters
    ----------
    receptor : str
        Key of the gene containing the molecule acting as receptor (protein)
    ligand : str
        Key of the gene containing the molecule acting as ligand

    Returns
    -------
    float
        Interaction energy in kcal/mol, as reported by AutoDock Vina --score-only.
    """
    _validate = {
        parse.Required('receptor'): parse.Molecule_name,
        parse.Required('ligand'): parse.Molecule_name,
        }

    def __init__(self, receptor='Protein', ligand='Ligand', *args, **kwargs):
        ObjectiveProvider.__init__(self, **kwargs)
        self.receptor = receptor
        self.ligand = ligand
        self._paths = []
        if os.name == 'posix' and os.path.exists('/dev/shm'):
            self.tmpdir = '/dev/shm'
        else:
            self.tmpdir = _get_default_tempdir()

    def evaluate(self, ind):
        """
        Run a subprocess calling Vina binary with provided options,
        and parse the results. Clean tmp files at exit.
        """
        self.tmpfile = os.path.join(self.tmpdir, next(_get_candidate_names()))
        receptor = ind.find_molecule(self.receptor)
        ligand = ind.find_molecule(self.ligand)

        receptor_pdbqt = self.prepare_receptor(receptor)
        ligand_pdbqt = self.prepare_ligand(ligand)
        command = ['vina', '--score_only', '--cpu', '1',
                   '--receptor', receptor_pdbqt, '--ligand', ligand_pdbqt]

        try:
            stream = check_output(command, universal_newlines=True)
        except CalledProcessError:
            logger.warning("Could not run Vina with command %s", command)
            return -100000 * self.weight
        else:
            return self.parse_output(stream)
        finally:
            self.clean()

    def prepare_molecule(self, molecule, which='receptor'):
        """ Write PDBQT for molecule """
        if which == 'ligand':
            script = [sys.executable, prepare_ligand4.__file__.replace('.pyc', '.py'), '-l']
        else:
            script = [sys.executable, prepare_receptor4.__file__.replace('.pyc', '.py'), '-r']
        print(script)
        path = '{}_{}.pdb'.format(self.tmpfile, which)
        pathqt = path + 'qt'
        pdb = molecule.write(absolute=path, filetype='pdb')
        self._paths.append(path)
        output = call(script + [path, '-o', pathqt])
        self._paths.append(pathqt)
        return pathqt

    def prepare_receptor(self, molecule):
        path = '{}_receptor.pdb'.format(self.tmpfile)
        pathqt = path + 'qt'
        pdb = molecule.write(absolute=path, filetype='pdb')
        self._paths.append(path)
        mol = MolKit.Read(path)[0]
        mol.buildBondsByDistance()
        RPO = AD4ReceptorPreparation(mol, outputfilename=pathqt)
        self._paths.append(pathqt)
        return pathqt

    def prepare_ligand(self, molecule):
        path = '{}_ligand.pdb'.format(self.tmpfile)
        pathqt = path + 'qt'
        pdb = molecule.write(absolute=path, filetype='pdb')
        self._paths.append(path)
        mol = MolKit.Read(path)[0]
        mol.buildBondsByDistance()
        RPO = AD4LigandPreparation(mol, outputfilename=pathqt)
        #    inactivate_all_torsions=True)
        self._paths.append(pathqt)
        return pathqt

    def parse_output(self, stream):
        for line in stream.splitlines():
            if line[:9] == "Affinity:":
                return float(line.split()[1])
        return -1000 * self.weight

    def clean(self):
        for p in self._paths:
            os.remove(p)
        self._paths = []