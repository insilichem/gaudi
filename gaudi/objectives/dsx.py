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
This objective is a wrapper around the binaries provided by
Neudert and Klebe at http://pc1664.pharmazie.uni-marburg.de/drugscore/
and calculates the score of the current pose.

The lower, the better, so usually you will use a -1.0 weight.

"""

# Python
import os
import subprocess
import tempfile
import logging
# GAUDI
from gaudi.objectives import ObjectiveProvider
from gaudi import parse

logger = logging.getLogger(__name__)


def enable(**kwargs):
    return DSX(**kwargs)


class DSX(ObjectiveProvider):

    """
    DSX class

    Parameters
    ----------
    binary : str
        Path to the DSX binary
    potentials : str
        Path to DSX potentials
    protein : str
        The molecule name that is acting as a protein
    ligand : str
        The molecule name that is acting as a ligand
    terms : list of bool
        Enable or disable certain terms in the score function in
        this order:

            distance-dependent pair potentials, 
            torsion potentials,
            intramolecular clashes, 
            sas potentials, 
            hbond potentials

    sorting : int
        Sorting mode. An int between 0-6, read binary help for -S
    cofactor_handling : int
        Cofactor handling mode. An int between 0-7, read binary help for -I
    """
    validate = parse.Schema({
        parse.Required('binary'): parse.Path(isfile=True),
        parse.Required('potentials'): parse.Path(isdir=True),
        parse.Required('protein'): parse.Molecule_name,
        parse.Required('ligand'): parse.Molecule_name,
        parse.Required('terms'): parse.All([parse.Boolean], parse.Length(min=5, max=5)),
        'sorting': parse.All(parse.Coerce(int), parse.Range(min=0, max=6)),
        'cofactor_handling': parse.All(parse.Coerce(int), parse.Range(min=0, max=7))
        })
    def __init__(self, binary=None, potentials=None, protein='Protein',
                 ligand='Ligand', terms=None, sorting=1, cofactor_handling=1,
                 *args, **kwargs):
        ObjectiveProvider.__init__(self, **kwargs)
        self.binary = binary
        self.potentials = potentials
        self.protein_names = protein
        self.ligand_names = ligand
        self.terms = terms
        self.sorting = sorting
        self.cofactor_handling = cofactor_handling
        self.oldworkingdir = os.getcwd()
        self.tempdir = tempfile._get_default_tempdir()

    def get_molecule_by_name(self, ind, *names):
        """
        Get a molecule gene instance of individual by its name
        """
        for name in names:
            for gene in ind.genes.values():
                if gene.__class__.__name__ == 'Molecule' and gene.name == name:
                    yield gene

    def evaluate(self, ind):
        """
        Run a subprocess calling DSX binary with provided options,
        and parse the results. Clean tmp files at exit.
        """
        tmpfile = os.path.join(self.tempdir,
                               next(tempfile._get_candidate_names()))

        # Retrieve proteins and write them to single mol2 file
        proteinpath = '{}_protein.mol2'.format(tmpfile)
        proteins = list(self.get_molecule_by_name(ind, *self.protein_names))
        last_protein = proteins.pop()
        last_protein.write(absolute=proteinpath, combined_with=proteins)
        # Retrieve ligands and write them to single mol2 file
        ligandpath = '{}_ligand.mol2'.format(tmpfile)
        ligands = list(self.get_molecule_by_name(ind, *self.ligand_names))
        last_ligand = ligands.pop()
        last_ligand.write(absolute=ligandpath, combined_with=ligands)

        T0, T1, T2, T3, T4 = [1.0 * t for t in self.terms]
        command = map(str, (self.binary, '-P', proteinpath, '-L', ligandpath,
                            '-I', self.cofactor_handling, '-S', self.sorting,
                            '-T0', T0, '-T1', T1, '-T2', T2, '-T3', T3, '-T4',
                            T4, '-D', self.potentials))
        os.chdir(tempfile._get_default_tempdir())
        try:
            stream = subprocess.check_output(command, universal_newlines=True)
        except subprocess.CalledProcessError:
            logger.warning("Could not run DSX with command %s", command)
            return -100000 * self.weight
        else:
            # 1. Get output filename from stdout (located at working directory)
            # 2. Find line '@RESULTS' and go to sixth line below
            # 3. The score is in the first row of the table, at the third field
            dsx_results = os.path.join(os.getcwd(),
                                       stream.splitlines()[-2].split()[-1])
            with open(dsx_results) as f:
                lines = f.read().splitlines()
                i = lines.index('@RESULTS')
                score = lines[i + 4].split('|')[3].strip()
                return float(score)
        finally:
            os.remove(dsx_results)
            os.remove(proteinpath)
            os.remove(ligandpath)
            os.chdir(self.oldworkingdir)
