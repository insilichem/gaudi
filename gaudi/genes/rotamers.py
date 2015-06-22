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
:mod:`gaudi.genes.rotamers` allows to explore side chains flexibility
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
import os
# Chimera
from Rotamers import getRotamers, useRotamer, NoResidueRotamersError
import SwapRes
# External dependencies
from repoze.lru import LRUCache
import deap.tools
# GAUDI
from gaudi.genes import GeneProvider
from gaudi.parse import parse_rawstring

logger = logging.getLogger(__name__)


class Rotamers(GeneProvider):

    def __init__(self, residues=None, library='Dunbrack',
                 mutations=[], ligation=False,
                 **kwargs):
        GeneProvider.__init__(self, **kwargs)
        self._kwargs = kwargs
        self._residues = residues
        self.library = library
        self.mutations = mutations
        self.ligation = ligation
        self.allele = []
        # set caches
        if self.name + '_res' not in self._cache:
            self._cache[self.name + '_res'] = OrderedDict()
        if self.name + '_rot' not in self._cache:
            self._cache[self.name + '_rot'] = LRUCache(300)

        if self.ligation:
            self.random_number = random.random()
        else:
            self.random_number = None

    @property
    def residues(self):
        return self._cache[self.name + '_res']

    @property
    def rotamers(self):
        return self._cache[self.name + '_rot']

    def __ready__(self):
        self._residues_rawstring = tuple(
            parse_rawstring(r) for r in self._residues)
        for molecule, resid in self._residues_rawstring:
            try:
                res = next(r for r in self.parent.genes[molecule].compound.mol.residues
                           if r.id.position == resid)
            # molecule or residue not found
            except (KeyError, StopIteration):
                raise
            else:  # residue was found!
                self.residues[(molecule, resid)] = res
                self.allele.append((self.choice(self.mutations + [res.type]),
                                    random.random()))

    def express(self):
        for (mol, pos), (restype, i) in zip(self.residues, self.allele):
            try:
                rot = self.get_rotamers(mol, pos, restype)
            except NoResidueRotamersError:  # ALA, GLY...
                SwapRes.swap(self.residues[(mol, pos)], restype)
            else:
                useRotamer(self.residues[(mol, pos)], [rot[int(i * len(rot))]])
            finally:
                self.residues[(mol, pos)] = res = \
                    next(r for r in self.parent.genes[mol].compound.mol.residues
                         if r.id.position == pos)
                for a in res.atoms:
                    a.display = 1

    def unexpress(self):
        for res in self.residues.values():
            for a in res.atoms:
                a.display = 0

    def mate(self, mate):
        if self.ligation:
            self_residues, self_rotamers = zip(*self.allele)
            mate_residues, mate_rotamers = zip(*mate.allele)
            self_rotamers, mate_rotamers = deap.tools.cxTwoPoint(
                list(self_rotamers), list(mate_rotamers))
            self.allele = map(list, zip(self_residues, self_rotamers))
            mate.allele = map(list, zip(mate_residues, mate_rotamers))
        else:
            self.allele, mate.allele = deap.tools.cxTwoPoint(
                self.allele, mate.allele)

    def mutate(self, indpb):
        if random.random() < self.indpb:
            self.allele[:] = []
            if self.ligation:  # don't forget to get a new random!
                self.random_number = random.random()
            for res in self.residues.values():
                self.allele.append(
                    (self.choice(self.mutations + [res.type]),
                        random.random()
                     )
                )

    ###

    def choice(self, l):
        """
        Overrides ``random.choice`` with custom one so we can
        reuse a previously obtained random number. This helps dealing
        with the ``ligation`` parameter, which forces all the requested
        residues to mutate to the same type
        """
        if self.random_number:
            return l[int(self.random_number * len(l))]
        return l[int(random.random() * len(l))]

    def get_rotamers(self, mol, pos, restype):
        rotamers = self.rotamers.get((mol, pos, restype))
        if not rotamers:
            try:
                rotamers = getRotamers(self.residues[(mol, pos)], resType=restype,
                                       lib=self.library.title())[1]
            except NoResidueRotamersError:  # ALA, GLY... has no rotamers
                raise
            except KeyError:
                raise
            else:
                self.rotamers.put((mol, pos, restype), rotamers)
        return rotamers


def enable(**kwargs):
    return Rotamers(**kwargs)
