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
# Chimera
import chimera
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
        parse.Required('binary'): parse.ExpandUserPathExists,
        parse.Required('potentials'): parse.ExpandUserPathExists,
        parse.Required('proteins'): [parse.Molecule_name],
        parse.Required('ligands'): [parse.Molecule_name],
        parse.Required('terms'): parse.All([parse.Boolean()], parse.Length(min=5, max=5)),
        'sorting': parse.All(parse.Coerce(int), parse.Range(min=0, max=6)),
        'cofactor_mode': parse.All(parse.Coerce(int), parse.Range(min=0, max=7)),
        'with_covalent': parse.Coerce(bool),
        'with_metals': parse.Coerce(bool)
        }, extra=parse.ALLOW_EXTRA)
    
    def __init__(self, binary=None, potentials=None, proteins=('Protein',),
                 ligands=('Ligand',), terms=None, sorting=1, cofactor_mode=0,
                 with_covalent=False, with_metals=True, *args, **kwargs):
        ObjectiveProvider.__init__(self, **kwargs)
        self.binary = binary
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
            self.tmpdir = tempfile._get_default_tempdir()

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
        self.tmpfile = os.path.join(self.tmpdir, next(tempfile._get_candidate_names()))
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
        # Retrieve proteins and write them to single mol2 file
        proteinpath = '{}_proteins.pdb'.format(self.tmpfile)
        last_protein = proteins.pop()
        last_protein.write_pdb(absolute=proteinpath, combined_with=proteins)
        self._paths['proteins'] = proteinpath

    def prepare_ligands(self, ligands):
        # Retrieve ligands and write them to single mol2 file
        ligandpath = '{}_ligands.mol2'.format(self.tmpfile)
        metalpath = '{}_metals.mol2'.format(self.tmpfile)
        ligand_mols = [lig.compound.mol for lig in ligands]

        if self.with_metals:
            # Split metals from ligand
            nonmetal_mols, metal_mols = [], []
            for ligand in ligands:
                nonmetals, metals = [], []
                for atom in ligand.compound.mol.atoms:
                    if atom.element in chimera.elements.metals:
                        metals.append(atom)
                    else:
                        nonmetals.append(atom)
                nonmetal_mols.append(molecule_from_atoms(ligand.compound.mol, nonmetals))
                if metals:
                    metal_mols.append(molecule_from_atoms(ligand.compound.mol, metals))
            if metal_mols:
                writeMol2(metal_mols, metalpath, temporary=True)
                self._paths['metals'] = metalpath
                ligand_mols = nonmetal_mols
       
        writeMol2(ligand_mols, ligandpath, temporary=True)
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
        cmd.extend(['-I', self.cofactor_mode])
        cmd.extend(['-S', self.sorting])
        T0, T1, T2, T3, T4 = [1.0 * t for t in self.terms]
        cmd.extend(['-T0', T0, '-T1', T1, '-T2', T2, '-T3', T3, '-T4', T4])
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