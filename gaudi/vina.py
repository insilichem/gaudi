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
    prepare_each : bool
        Whether to prepare receptors and ligands in every evaluation or try
        to cache the results for faster performance.

    Returns
    -------
    float
        Interaction energy in kcal/mol, as reported by AutoDock Vina --score-only.

    Notes
    -----
    - AutoDock scripts ``prepare_ligand4.py`` and ``prepare_receptor4.py`` are
    used to prepare the corresponding .pdqt files that will be used as input for 
    AutoDock Vina scorer.
    - No repairs nor cleanups will be performed on ligand/receptor molecules, so 
    the user has to take into account that provided .mol2 or .pdb files have 
    correct atom types and correct structure (including Hydrogen atoms that will 
    be taken into account in the docking evaluation). Otherwise, AutoDock 
    errors/warnings could appear (e.g. ``ValueError: Could not find atomic number 
    for Lp Lp``)
    - Gasteiger charges will be added during the preparation of the .pdbqt files.
    - All torsions of the ligand will be marked as ``inactive`` for AutoDock, 
    because torsion changes are part of GaudiMM genes.
    """
    _validate = {
        parse.Required('receptor'): parse.Molecule_name,
        parse.Required('ligand'): parse.Molecule_name,
        'prepare_each': bool,
        }

    def __init__(self, receptor='Protein', ligand='Ligand', prepare_each=False,
                 *args, **kwargs):
        ObjectiveProvider.__init__(self, **kwargs)
        self.receptor = receptor
        self.ligand = ligand
        self.prepare_each = prepare_each
        self._paths = []
        self._tmpfile = None
        if os.name == 'posix' and os.path.exists('/dev/shm'):
            self.tmpdir = '/dev/shm'
        else:
            self.tmpdir = _get_default_tempdir()

    def evaluate(self, ind):
        """
        Run a subprocess calling Vina binary with provided options,
        and parse the results. Clean tmp files at exit.
        """
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

    def _prepare(self, molecule, which='receptor'):
        if which == 'receptor':
            preparer = AD4ReceptorPreparation
            kwargs = {"repairs": '', "cleanup": ''}
        elif which == 'ligand':
            preparer = AD4LigandPreparation
            kwargs = {"repairs": '', "cleanup": '', "inactivate_all_torsions": True}
        else:
            raise ValueError('which must be receptor or ligand')
        path = '{}_{}.pdb'.format(self.tmpfile, which)
        pathqt = path + 'qt'
        if not os.path.isfile(pathqt):
            pdb = molecule.write(absolute=path, filetype='pdb')
            self._paths.append(path)
            mol = MolKit.Read(path)[0]
            mol.buildBondsByDistance()
            RPO = preparer(mol, outputfilename=pathqt, **kwargs)
            self._paths.append(pathqt)
        else:
            # update coordinates
            self._update_pdbqt_coordinates(molecule.xyz(), pathqt)
        return pathqt

    def prepare_receptor(self, molecule):
        return self._prepare(molecule, 'receptor')

    def prepare_ligand(self, molecule):
        return self._prepare(molecule, 'ligand')

    @staticmethod
    def _update_pdbqt_coordinates(xyz, path):
        def is_atom(line):
            if line[0:6] in ('ATOM  ', 'HETATM'):
                for n in (line[30:38], line[38:46], line[46:54]):
                    try:
                        float(n)
                    except:
                        return False
                return True
            return False

        with open(path, 'r+') as f:
            lines = []
            i = 0
            for line in f:
                if is_atom(line):
                    line = line[:30] + '{:8.3f}{:8.3f}{:8.3f}'.format(*xyz[i]) + line[54:]
                    i += 1
                lines.append(line)
            f.seek(0)
            f.write(''.join(lines))
            f.truncate()


    def parse_output(self, stream):
        for line in stream.splitlines():
            if line[:9] == "Affinity:":
                return float(line.split()[1])
        return -1000 * self.weight

    def clean(self):
        if not self.prepare_each:
            return
        for p in self._paths:
            os.remove(p)
        self._paths = []

    @property
    def tmpfile(self):
        if self.prepare_each or self._tmpfile is None:
            self._tmpfile = os.path.join(self.tmpdir, next(_get_candidate_names()))
        return self._tmpfile
