#!/usr/bin/env python
# -*- coding: utf-8 -*-

##############
# GAUDIasm: Genetic Algorithms for Universal
# Design Inference and Atomic Scale Modeling
# Authors:  Jaime Rodriguez-Guerra Pedregal
#            <jaime.rodriguezguerra@uab.cat>
#           Jean-Didier Marechal
#            <jeandidier.marechal@uab.cat>
# Web: https://bitbucket.org/insilichem/gaudi
##############

"""
This modules allows to explore side chains flexibility
in proteins, as well as mutation.

It needs that at least a :class:`gaudi.genes.rotamers.molecule.Molecule` has been
requested in the input file. Residues of those are referenced in the `residues` argument.

It also allows mutations in the selected residues. However, the resulting structure keeps
the same backbone, which may not be representative of the in-vivo behaviour. Use with caution.

"""

# Python
import random
from collections import OrderedDict
import logging
# Chimera
import chimera
from chimera import BondRot, dihedral
from chimera.phipsi import chiAtoms, AtomsMissingError
from Rotamers import getRotamerParams, NoResidueRotamersError
# External dependencies
import deap.tools
# GAUDI
from gaudi import parse
from gaudi.genes import GeneProvider

logger = logging.getLogger(__name__)


def enable(**kwargs):
    kwargs = Rotamers.validate(kwargs)
    return Rotamers(**kwargs)


class Rotamers(GeneProvider):

    """
    Rotamers class

    Parameters
    ----------
    residues : list of str
        Residues that should be analyzed. This has to be in the form:

            [ Protein/233, Protein/109 ]

        where the first element (before slash) is the gaudi.genes.molecule name
        and the second element (after slash) is the residue position number in that
        molecule.

        This list of str is later parsed to the proper chimera.Residue objects

    library : {'Dunbrack', 'Dynameomics'}
        The rotamer library to use.
    """
    _validate = {
        parse.Required('residues'): [parse.Named_spec("molecule", "residue")],
        'library': parse.Any('Dunbrack', 'dunbrack', 'Dynameomics', 'dynameomics'),
       }

    # Avoid unnecesary calls to expensive get_rotamers if residue is known
    # to not have any rotamers
    _residues_without_rotamers = set(('ALA', 'GLY'))

    def __init__(self, residues=None, library='Dunbrack', **kwargs):
        GeneProvider.__init__(self, **kwargs)
        self._kwargs = kwargs
        self._residues = residues
        self.library = library
        self.allele = []
        # set caches
        self.residues = self._cache.setdefault(self.name + '_residues', OrderedDict())

    def __ready__(self):
        """
        Second stage of initialization.

        It parses the requested residues strings to actual residues.
        """
        for molname, pos in self._residues:
            residues = self.parent.find_molecule(molname).find_residues(pos)
            if len(residues) > 1 and pos != '*':
                logger.warn('Found one more than residue for %s/%s', molname, pos)
            for r in residues:
                self.patch_residue(r)
                self.residues[(molname, r.id.position)] = r
                self.allele.append(random.random())
    
    def __deepcopy__(self, memo):
        new = self.__class__(residues=self._residues, library=self.library, **self._kwargs)
        new.residues = self.residues
        new.allele = self.allele[:]
        return new

    def express(self):
        for ((molname, pos), residue), i in zip(self.residues.items(), self.allele):
            if residue.type not in self._residues_without_rotamers:
                try:
                    rotamers = getRotamerParams(residue, lib=self.library.title())[2]
                except NoResidueRotamersError:  # ALA, GLY...
                    logger.warn('%s/%s (%s) has no rotamers', molname, pos, residue.type)
                    self._residues_without_rotamers.add(residue.type)
                else:
                    rotamer = rotamers[int(i * len(rotamers))]
                    self.update_rotamer(residue, rotamer.chis)

    def unexpress(self):
        for res in self.residues.values():
            for torsion in res._rotamer_torsions:
                torsion.reset()

    def mate(self, mate):
        self.allele, mate.allele = deap.tools.cxTwoPoint(self.allele, mate.allele)

    def mutate(self, indpb):
        self.allele = [random.random() if random.random() < indpb else i for i in self.allele]

    @staticmethod
    def update_rotamer(residue, chis):
        for bondrot, chi in zip(residue._rotamer_torsions, chis):
            bondrot.adjustAngle(chi - bondrot.chi, bondrot.biggerSide())

    @staticmethod
    def patch_residue(residue):
        if getattr(residue, '_rotamer_torsions', None):
            return
        residue._rotamer_torsions = [] # BondRot objects cache
        for chi in range(1, 5):
            try:
                atoms = chiAtoms(residue, chi)
                bond = atoms[1].bondsMap[atoms[2]]
                br = BondRot(bond)
                br.anchor = br.biggerSide()
                br.chi = dihedral(*[a.coord() for a in atoms])
                residue._rotamer_torsions.append(br)
            except AtomsMissingError:
                break
            except (chimera.error, ValueError) as v:
                    if "cycle" in str(v) or "already used" in str(v):
                        continue  # discard bonds in cycles and used!
                    break
    @staticmethod
    def all_chis(residue):
        chis = []
        for i in range(1, 5):
            try:
                chis.append(getattr(residue, 'chi{}'.format(i)))
            except AttributeError:
                break
        return chis