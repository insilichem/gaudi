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
    _residues_without_rotamers = ['ALA', 'GLY']

    def __init__(self, residues=None, library='Dunbrack', **kwargs):
        GeneProvider.__init__(self, **kwargs)
        self._kwargs = kwargs
        self._residues = residues
        self.library = library
        self.allele = []
        # set caches
        self.residues = self._cache.setdefault(self.name + '_residues', OrderedDict())
        self.rotamers = self._cache.setdefault(self.name + 'rotamers', dict())

    def __ready__(self):
        """
        Second stage of initialization.

        It parses the requested residues strings to actual residues.
        """
        for molname, pos in self._residues:
            try:
                mol = self.parent._molecules[molname].compound.mol
            except KeyError:
                available_molnames = self.parent._molecules.keys()
                raise KeyError('Molecule {} was not found!'
                               'Try with one of {}'.format(molname, available_molnames))
            if pos == '*':
                residues = mol.residues
            else:
                residues = [r for r in mol.residues if r.id.position == pos]
                if len(residues) > 1:
                    logger.warn('Found one more than residue for %s/%s', molname, pos)
            for r in residues:
                r._original_chis = [] # Save torsion angles of initial rotamer
                for i in range(1, 5):
                    try:
                        chi = getattr(r, 'chi' + str(i))
                        r._original_chis.append(chi)
                    except AttributeError:
                        continue 
                self.residues[(molname, r.id.position)] = r
                self.allele.append(random.random())
    
    def __deepcopy__(self, memo):
        new = self.__class__(residues=self._residues, library=self.library, **self._kwargs )
        new.residues = self.residues
        new.rotamers = self.rotamers
        new.allele = self.allele[:]
        return new

    def express(self):
        for ((molname, pos), residue), i in zip(self.residues.items(), self.allele):
            try:
                rotamers = self.get_rotamers((molname, pos), residue)
            except NoResidueRotamersError:  # ALA, GLY...
                logger.warn('%s/%s (%s) has no rotamers', molname, pos, residue.type)
            else:
                rotamer = rotamers[int(i * len(rotamers))]
                self.update_rotamer(residue, rotamer.chis)

    def unexpress(self):
        for res in self.residues.values():
            self.update_rotamer(res, res._original_chis)

    def mate(self, mate):
        self.allele, mate.allele = deap.tools.cxTwoPoint(self.allele, mate.allele)

    def mutate(self, indpb):
        self.allele = [random.random() if random.random() < indpb else i for i in self.allele]

    def get_rotamers(self, key, residue):
        """
        Gets the requested rotamers out of cache and if not found,
        creates the library of chis and stores it in the cache.

        Parameters
        ----------
        residue : chimera.Residue
            The residue that must be analyzed

        Returns
        -------
            List of Rotamer objects
        """
        rotamers = getRotamerParams(residue, lib=self.library.title())[2]
        return self.rotamers.setdefault(key, rotamers)

    @staticmethod
    def update_rotamer(residue, chis):
        for i, chi in enumerate(chis):
            setattr(residue, 'chi{}'.format(i + 1), chi)
