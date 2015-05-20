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

    def __init__(self, residues=None, library='Dunbrack', mutations=[],
                 **kwargs):
        GeneProvider.__init__(self, **kwargs)
        self._kwargs = kwargs
        self._residues = residues
        self.library = library
        self.mutations = mutations
        self.allele = []
        # set (or retrieve) caches
        try:
            self.residues = self._cache['residues']
            self.rotamers = self._cache['rotamers']
        except KeyError:
            self.residues = self._cache['residues'] = OrderedDict()
            self.rotamers = self._cache['rotamers'] = LRUCache(300)

        # find requested residues
        self._residues_rawstring = tuple(parse_rawstring(r) for r in residues)
        for molecule, resid in self._residues_rawstring:
            try:
                res = next(r for r in self.parent.genes[molecule].compound.mol.residues
                           if r.id.position == resid)
            except (KeyError, StopIteration):  # molecule or residue not found
                raise
            else:  # residue was found!
                self.residues[(molecule, resid)] = res
                self.allele.append(
                    (random.choice(self.mutations + [res.type]),
                        random.random()
                     )
                )

    def __deepcopy__(self, memo):
        new = self.__class__(self._residues, self.library, self.mutations,
                             **self._kwargs)
        new.__dict__.update((k, v) for k, v in self.__dict__.items())
        new.allele = self.allele[:]
        return new

    def express(self):
        for (mol, pos), (restype, i) in zip(self.residues, self.allele):
            try:
                rot = self.get_rotamers(mol, pos, restype)
            except NoResidueRotamersError:  # ALA, GLY...
                SwapRes.swap(self.residues[(mol, pos)], restype, preserve=True)
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
        self.allele, mate.allele = deap.tools.cxTwoPoint(
            self.allele, mate.allele)

    def mutate(self, indpb):
        if random.random() < self.indpb:
            self.allele = []
            for molecule, resid in self._residues_rawstring:
                try:
                    res = next(r for r in self.parent.genes[molecule].compound.mol.residues
                               if r.id.position == resid)
                # molecule or residue not found
                except (KeyError, StopIteration):
                    raise
                else:  # residue was found!
                    self.residues[(molecule, resid)] = res
                    self.allele.append(
                        (random.choice(self.mutations + [res.type]),
                            random.random()
                         )
                    )

    def write(self, path, name):
        pass
        # rotamerline = '{}.{} {} {} {}\n'.format(res.id.position, res.id.chainId,
        # self.library.title(), res.type, ' '.join(map(str,rot.chis)))

    ###
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
