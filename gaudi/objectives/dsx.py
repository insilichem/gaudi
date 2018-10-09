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
This objective is a wrapper around the binaries `provided by
Neudert and Klebe <http://pc1664.pharmazie.uni-marburg.de/drugscore/>`_
and calculates the score of the current pose.

The lower, the better, so usually you will use a -1.0 weight.

"""

# Python
import os
import subprocess
from tempfile import _get_default_tempdir as default_tempdir, _get_candidate_names as tempnames
from distutils.spawn import find_executable
import logging
# Chimera
from SplitMolecule.split import molecule_from_atoms
from WriteMol2 import writeMol2
# GAUDI
from gaudi.objectives import ObjectiveProvider
from gaudi import parse

logger = logging.getLogger(__name__)


def enable(**kwargs):
    kwargs = DSX.validate(kwargs)
    return DSX(**kwargs)


class DSX(ObjectiveProvider):

    """
    DSX class

    Parameters
    ----------
    protein : str
        The molecule name that is acting as a protein
    ligand : str
        The molecule name that is acting as a ligand
    binary : str, optional
        Path to the DSX binary. Only needed if ``drugscorex`` is not in PATH.
    potentials : str, optional
        Path to DSX potentials. Only needed if ``DSX_POTENTIALS`` env var has not
        been set by the installation process (``conda install -c insilichem drugscorex``
        normally takes care of that).
    terms : list of bool, optional
        Enable (True) or disable (False) certain terms in the score function in
        this order: distance-dependent pair potentials, torsion potentials,
        intramolecular clashes, sas potentials, hbond potentials

    sorting : int, defaults to 1
        Sorting mode. An int between 0-6, read binary help for -S::

            -S int :  Here you can specify the mode that affects how the results
                    will be sorted. The default mode is '-S 1', which sorts the
                    ligands in the same order as they are found in the lig_file.
                    The following modes are possible::

                        0: Same order as in the ligand file
                        1: Ordered by increasing total score
                        2: Ordered by increasing per-atom-score
                        3: Ordered by increasing per-contact-score
                        4: Ordered by increasing rmsd
                        5: Ordered by increasing torsion score
                        6: Ordered by increasing per-torsion-score

    cofactor_mode : int, defaults to 0
        Cofactor handling mode. An int between 0-7, read binary help for -I::

            -I int :  Here you can specify the mode that affects how cofactors,
                    waters and metals will be handeled.
                    The default mode is '-I 1', which means, that all molecules
                    are treated as part of the protein. If a structure should
                    not be treated as part of the protein you have supply a
                    seperate file with seperate MOLECULE entries corresponding
                    to each MOLECULE entry in the ligand_file (It is assumed
                    that the structure, e.g. a cofactor, was kept flexible in
                    docking, so that there should be a different geometry
                    corresponding to each solution. Otherwise it won't make
                    sense not to treat it as part of the protein.).
                    The following modes are possible:
                        0: cofactors, waters and metals interact with protein,
                        ligand and each other
                        1: cofactors, waters and metals are treated as part of
                        the protein
                        2: cofactors and metals are treated as part of the protein
                        (waters as in mode 0)
                        3: cofactors and waters are treated as part of the protein
                        4: cofactors are treated as part of the protein
                        5: metals and waters are treated as part of the protein
                        6: metals are treated as part of the protein
                        7: waters are treated as part of the protein
                    Please note: Only those structures can be treated
                    individually, which are supplied in seperate files.
    with_covalent : bool, defaults to False
        Whether to deal with covalently bonded atoms as normal atoms (False) or not (True)
    with_metals : bool, defaults to True
        Whether to deal with metal atoms as normal atoms (False) or not (True)

    Returns
    -------
    float
        Interaction energy as reported by DSX output logs.
    """
    _validate = {
        parse.Required('proteins'): [parse.Molecule_name],
        parse.Required('ligands'): [parse.Molecule_name],
        'binary': parse.ExpandUserPathExists,
        'potentials': parse.ExpandUserPathExists,
        'terms': parse.All([parse.Coerce(bool)], parse.Length(min=5, max=5)),
        'sorting': parse.All(parse.Coerce(int), parse.Range(min=0, max=6)),
        'cofactor_mode': parse.All(parse.Coerce(int), parse.Range(min=0, max=7)),
        'with_covalent': parse.Coerce(bool),
        'with_metals': parse.Coerce(bool)
        }

    def __init__(self, binary=None, potentials=None, proteins=('Protein',),
                 ligands=('Ligand',), terms=None, sorting=1, cofactor_mode=0,
                 with_covalent=False, with_metals=True, *args, **kwargs):
        ObjectiveProvider.__init__(self, **kwargs)
        self.binary = find_executable('drugscorex') if binary is None else binary
        if not self.binary:
            raise ValueError('Could not find `drugscorex` executable. Please install it '
                             'with `conda install -c insilichem drugscorex` or manually '
                             'specify the location with `binary` and `potentials` keys.')
        self.potentials = potentials
        self.protein_names = proteins
        self.ligand_names = ligands
        self.terms = terms
        self.sorting = sorting
        self.cofactor_mode = cofactor_mode
        self.with_covalent = with_covalent
        self.with_metals = with_metals

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
        self.tmpfile = os.path.join(self.tmpdir, next(tempnames()))
        proteins = list(self.get_molecule_by_name(ind, *self.protein_names))
        ligands = list(self.get_molecule_by_name(ind, *self.ligand_names))

        self.prepare_proteins(proteins)
        self.prepare_ligands(ligands)
        command = self.prepare_command()

        try:
            os.chdir(self.tmpdir)
            stream = subprocess.check_output(command, universal_newlines=True)
        except subprocess.CalledProcessError:
            logger.warning("Could not run DSX with command %s", command)
            return -100000 * self.weight
        else:
            return self.parse_output(stream)
        finally:
            self.clean()
            os.chdir(self._oldworkingdir)

    def prepare_proteins(self, proteins):
        proteinpath = '{}_proteins.pdb'.format(self.tmpfile)
        last_protein = proteins.pop()
        last_protein.write(absolute=proteinpath, combined_with=proteins, filetype='pdb')
        self._paths['proteins'] = proteinpath

    def prepare_ligands(self, ligands):
        ligandpath = '{}_ligands.mol2'.format(self.tmpfile)
        metalpath = '{}_metals.mol2'.format(self.tmpfile)
        ligand_mols = [lig.compound.mol for lig in ligands]

        if self.with_metals:
            # Split metals from ligand
            nonmetal_mols, metal_mols = [], []
            for ligand in ligand_mols:
                nonmetals, metals = [], []
                for atom in ligand.atoms:
                    if atom.element.isMetal:
                        metals.append(atom)
                    else:
                        nonmetals.append(atom)
                nonmetal_mols.append(molecule_from_atoms(ligand, nonmetals))
                if metals:
                    metal_mols.append(molecule_from_atoms(ligand, metals))
            if metal_mols:
                writeMol2(metal_mols, metalpath, temporary=True)
                self._paths['metals'] = metalpath
                ligand_mols = nonmetal_mols

        writeMol2(ligand_mols, ligandpath, temporary=True, multimodelHandling='combined')
        self._paths['ligands'] = ligandpath

    def prepare_command(self):
        cmd = [self.binary]
        if self.with_covalent:
            cmd.append('-c')
        cmd.extend(['-P', self._paths['proteins'], '-L', self._paths['ligands']])
        if self.with_metals:
            metalpath = self._paths.get('metals')
            if metalpath:
                cmd.extend(['-M', metalpath])
        if self.cofactor_mode is not None:
            cmd.extend(['-I', self.cofactor_mode])
        if self.sorting is not None:
            cmd.extend(['-S', self.sorting])
        if self.terms is not None:
            T0, T1, T2, T3, T4 = [1.0 * t for t in self.terms]
            cmd.extend(['-T0', T0, '-T1', T1, '-T2', T2, '-T3', T3, '-T4', T4])
        if self.potentials is not None:
            cmd.extend(['-D', self.potentials])
        return map(str, cmd)

    def parse_output(self, stream):
        # 1. Get output filename from stdout (located at working directory)
        # 2. Find line '@RESULTS' and go to sixth line below
        # 3. The score is in the first row of the table, at the third field
        dsx_results = os.path.join(self.tmpdir, stream.splitlines()[-2].split()[-1])
        self._paths['output'] = dsx_results
        with open(dsx_results) as f:
            lines = f.read().splitlines()
            i = lines.index('@RESULTS')
            score = lines[i + 4].split('|')[3].strip()
            return float(score)

    def clean(self):
        for p in self._paths.values():
            os.remove(p)
        self._paths.clear()
