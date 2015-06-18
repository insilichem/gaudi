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
:mod:`gaudi.genes.torsion` helps explore small molecules flexibility.

It does so by performing bond rotations in the selected :class:`gaudi.genes.molecule.Molecule`
objects.
"""

# Python
import random
import logging
from itertools import izip
from collections import OrderedDict
# Chimera
import chimera
# External dependencies
from deap.tools import cxSimulatedBinaryBounded, mutPolynomialBounded
# GAUDI
from gaudi.genes import GeneProvider
from gaudi import box

logger = logging.getLogger(__name__)


class Torsion(GeneProvider):
    BONDS_ROTS = {}

    def __init__(self, target=None, flexibility=None, max_bonds=30, **kwargs):
        GeneProvider.__init__(self, **kwargs)
        self._kwargs = kwargs
        self.target = target
        self.flexibility = 360.0 if flexibility > 360 else flexibility
        self.max_bonds = max_bonds
        self.nonrotatable = ()
        self.allele = [self.random_angle() for i in xrange(self.max_bonds)]

    def __ready__(self):
        if self.target not in self.parent.genes:
            logger.error("Gene for target %s is not in individual",
                         self.target)
            raise

        self.anchor = self._get_anchor()
        self.rotatable_bonds = list(self.get_rotatable_bonds())

    def express(self):
        self.rotatable_bonds[:] = []
        for alpha, br in izip(self.allele, self.get_rotatable_bonds()):
            try:
                if all(a.idatmType in ('C2', 'N2') for a in br.bond.atoms):
                    alpha = 0 if alpha < 180 else 180
                br.adjustAngle(alpha - br.angle, br.rotanchor)
            # A null bondrot was returned -> non-rotatable bond
            except AttributeError:
                continue
            else:
                self.rotatable_bonds.append(br)

    def unexpress(self):
        for br in self.rotatable_bonds:
            br.adjustAngle(-br.angle, br.rotanchor)

    def mate(self, mate):
        self.allele[:], mate.allele[:] = cxSimulatedBinaryBounded(
            self.allele, mate.allele, eta=self.cxeta,
            low=-0.5 * self.flexibility, up=0.5 * self.flexibility)

    def mutate(self, indpb):
        self.allele, = mutPolynomialBounded(self.allele,
                                            indpb=self.indpb, eta=self.mteta,
                                            low=-0.5 *
                                            self.flexibility,
                                            up=0.5 * self.flexibility)

    #####
    def random_angle(self):
        return random.uniform(-0.5 * self.flexibility, 0.5 * self.flexibility)

    def get_rotatable_bonds(self):
        atoms = self.parent.genes[self.target].compound.mol.atoms
        bonds = set(b for a in atoms for b in a.bonds if not a.element.isMetal)
        bonds = sorted(
            bonds, key=lambda b: min(y.serialNumber for y in b.atoms))

        self.anchor = self._get_anchor()
        for b in bonds:
            try:
                br = self.BONDS_ROTS[b]
            except KeyError:
                a1, a2 = b.atoms
                if a1 not in self.nonrotatable and \
                        a1.idatmType in ('C3', 'N3', 'C2', 'N2') and \
                        (a1.numBonds > 1 and a2.numBonds > 1) or \
                        a1.name == 'DUM' or a2.name == 'DUM':
                    try:
                        br = chimera.BondRot(b)
                    except (chimera.error, ValueError), v:
                        if "cycle" in str(v):
                            continue  # discard bonds in cycles!
                        elif "already used" in str(v):
                            print str(v)
                            continue
                        else:
                            raise
                    else:
                        br.rotanchor = box.find_nearest(self.anchor, b.atoms)
                        self.BONDS_ROTS[b] = br
                else:
                    continue
            yield br

    def update_rotatable_bonds(self):
        self.rotatable_bonds[:] = list(self.get_rotatable_bonds())

    def _get_anchor(self):
        try:
            search = next(g for g in self.parent.genes.values()
                          if g.__class__.__name__ == 'Search'
                          and g.target == self.target)
        except StopIteration:
            anchor = self.parent.genes[self.target].compound.donor
        else:
            try:
                anchor = next(a for a in self.parent.genes[self.target].atoms
                              if a.serialNumber == search.anchor)
            except (StopIteration, AttributeError):
                anchor = self.parent.genes[self.target].compound.donor

        anchor.name = 'ANC'
        return anchor

    def __deepcopy__(self, memo):
        new = self.__class__(self.target, self.flexibility, self.max_bonds,
                             **self._kwargs)
        new.__ready__()
        # new.__dict__.update((k, v) for k, v in self.__dict__.items())
        new.allele[:] = self.allele[:]
        return new


def enable(**kwargs):
    return Torsion(**kwargs)
