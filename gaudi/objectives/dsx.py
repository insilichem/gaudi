#!/usr/bin/python

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
:mod:`gaudi.objectives.angle` calculates the angle formed by three
given atoms (or the dihedral, if four atoms are given) and returns
the absolute difference of that angle and the target value.
"""

# Python
import os
import subprocess
import tempfile
import logging
# GAUDI
from gaudi.objectives import ObjectiveProvider

logger = logging.getLogger(__name__)


def enable(**kwargs):
    return DSX(**kwargs)


class DSX(ObjectiveProvider):

    def __init__(self, binary=None, potentials=None, protein='Protein',
                 ligand='Ligand', terms=None, sorting=1, cofactor_handling=1,
                 *args, **kwargs):
        ObjectiveProvider.__init__(self, **kwargs)
        self.binary = binary
        self.potentials = potentials
        self.protein_name = protein
        self.ligand_name = ligand
        self.terms = terms
        self.sorting = sorting
        self.cofactor_handling = cofactor_handling
        self.oldworkingdir = os.getcwd()

    @property
    def protein(self):
        for gene in self.parent.genes.values():
            if gene.__class__.__name__ == 'Molecule' and gene.name == self.protein_name:
                return gene

    @property
    def ligand(self):
        for gene in self.parent.genes.values():
            if gene.__class__.__name__ == 'Molecule' and gene.name == self.ligand_name:
                return gene

    def evaluate(self):
        tmpfile = os.path.join(tempfile._get_default_tempdir(),
                               next(tempfile._get_candidate_names()))
        proteinpath = '{}_protein.mol2'.format(tmpfile)
        ligandpath = '{}_ligand.mol2'.format(tmpfile)

        self.protein.write(absolute=proteinpath)
        self.ligand.write(absolute=ligandpath)

        T0, T1, T2, T3 = [1.0 * t for t in self.terms]
        command = map(str, (self.binary, '-P', proteinpath, '-L', ligandpath,
                            '-I', self.cofactor_handling, '-S', self.sorting,
                            '-T0', T0, '-T1', T1, '-T2', T2, '-T3', T3,
                            '-D', self.potentials))
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
            with open(os.path.join(os.getcwd(), stream.splitlines()[-2].split()[-1])) as f:
                lines = f.read().splitlines()
                i = lines.index('@RESULTS')
                score = lines[i + 4].split('|')[3].strip()
                return float(score)
        finally:
            os.chdir(self.oldworkingdir)
