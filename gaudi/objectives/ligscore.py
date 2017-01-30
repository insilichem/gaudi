#!/usr/bin/env python
# -*- coding: utf-8 -*-

##############
# GaudiMM: Genetic Algorithms with Unrestricted
# Descriptors for Intuitive Molecular Modeling
# Authors:  Jaime Rodriguez-Guerra Pedregal
#            <jaime.rodriguezguerra@uab.cat>
#           Jean-Didier Marechal
#            <jeandidier.marechal@uab.cat>
# Web: https://bitbucket.org/insilichem/gaudi
##############

"""
This objective is a wrapper around the binaries provided by
Neudert and Klebe at http://pc1664.pharmazie.uni-marburg.de/drugscore/
and calculates the score of the current pose.

The lower, the better, so usually you will use a -1.0 weight.

"""

# Python
from distutils.spawn import find_executable
import logging
import os
from tempfile import _get_default_tempdir as default_tempdir, _get_candidate_names as tempnames
import subprocess
# Chimera
from WriteMol2 import writeMol2
# GAUDI
from gaudi.objectives import ObjectiveProvider
from gaudi import parse

logger = logging.getLogger(__name__)


def enable(**kwargs):
    kwargs = LigScore.validate(kwargs)
    return LigScore(**kwargs)


class LigScore(ObjectiveProvider):

    """
    LigScore class

    Parameters
    ----------
    protein : list of str
        The name of molecules that are acting as proteins
    ligand : list of str
        The name of molecules that are acting as ligands
    binary : str, optional
        Path to ligand_score executable
    library : str, optional
        Path to LigScore lib file
    """
    _validate = {
        parse.Required('proteins'): [parse.Molecule_name],
        parse.Required('ligands'): [parse.Molecule_name],
        'method': parse.In(['rank', 'pose']),
        'binary': parse.ExpandUserPathExists,
        'library': parse.ExpandUserPathExists,
        }

    def __init__(self, proteins=('Protein',), ligands=('Ligand',), method='pose',
                 binary=None, library=None, *args, **kwargs):
        ObjectiveProvider.__init__(self, **kwargs)
        self.protein_names = proteins
        self.ligand_names = ligands
        self.binary = find_executable('ligand_score') if binary is None else binary
        self.library = library
        self.method = method

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
        Run a subprocess calling LigScore binary with provided options,
        and parse the results. Clean tmp files at exit.
        """
        self.tmpfile = os.path.join(self.tmpdir, next(tempnames()))
        proteins = list(self.get_molecule_by_name(ind, *self.protein_names))
        ligands = list(self.get_molecule_by_name(ind, *self.ligand_names))

        protein_path = self.prepare_proteins(proteins)
        ligand_path = self.prepare_ligands(ligands)
        command = self.prepare_command(protein_path, ligand_path)
        
        try:
            os.chdir(self.tmpdir)
            stream = subprocess.check_output(command, universal_newlines=True)
        except subprocess.CalledProcessError:
            logger.warning("Could not run LigScore with command %s", command)
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
        return proteinpath

    def prepare_ligands(self, ligands):
        ligandpath = '{}_ligands.mol2'.format(self.tmpfile)
        ligand_mols = [lig.compound.mol for lig in ligands]
       
        writeMol2(ligand_mols, ligandpath, temporary=True, multimodelHandling='combined')
        self._paths['ligands'] = ligandpath
        return ligandpath

    def prepare_command(self, protein_path, ligand_path):
        cmd = [self.binary, '--' + self.method, ligand_path, protein_path]
        if self.library:
            cmd.append(self.library)
        return map(str, cmd)

    def parse_output(self, stream):
        """ Get last word of first line (and unique) and parse it into float """
        return float(stream.splitlines()[0].split()[-1])

    def clean(self):
        for p in self._paths.values():
            os.remove(p)
        self._paths.clear()
