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
This module helps explore small molecules flexibility.

It does so by performing bond rotations in the selected `gaudi.genes.molecule.Molecule`
objects, if they exhibit free bond rotations.

"""

# Python
import random
import logging
from itertools import izip
# Chimera
import chimera
# External dependencies
from deap.tools import cxSimulatedBinaryBounded, mutPolynomialBounded
# GAUDI
from gaudi.genes import GeneProvider
from gaudi import box, parse

logger = logging.getLogger(__name__)


def enable(**kwargs):
    return Torsion(**kwargs)


class Torsion(GeneProvider):

    """
    Parameters
    ----------
    target: str
        Name of gaudi.genes.molecule instance to perform rotation on
    flexibility : int or float
        Maximum number of degrees a bond can rotate
    max_bonds :
        Expected number of free rotations in molecule. Needed to store
        arbitrary rotations. 

    Notes
    -----

    .. todo ::

        `max_bonds` should be automatically computed, based on ligand
        expected composition (careful with block-built ligands...)

    """

    validate = parse.Schema({
        parse.Required('target'): parse.Molecule_name,
        'flexibility': parse.Degrees,
        'max_bonds': parse.All(parse.Coerce(int), parse.Range(min=0))
        })

    BONDS_ROTS = {}

    def __init__(self, target=None, flexibility=None, max_bonds=30, **kwargs):
        GeneProvider.__init__(self, **kwargs)
        self._kwargs = kwargs
        self.target = target
        self.flexibility = 360.0 if flexibility > 360 else flexibility
        self.max_bonds = max_bonds
        self.nonrotatable = ()
        self.allele = [self.random_angle() for i in xrange(self.max_bonds)]

    def express(self):
        """
        Apply rotations to rotatable bonds
        """
        for alpha, br in izip(self.allele, self.rotatable_bonds):
            try:
                if all(a.idatmType in ('C2', 'N2') for a in br.bond.atoms):
                    alpha = 0 if alpha < 180 else 180
                br.adjustAngle(alpha - br.angle, br.rotanchor)
            # A null bondrot was returned -> non-rotatable bond
            except AttributeError:
                continue

    def unexpress(self):
        """
        Undo the rotations
        """
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
        """
        Returns a random angle within flexibility limits
        """
        return random.uniform(-0.5 * self.flexibility, 0.5 * self.flexibility)

    @property
    def rotatable_bonds(self):
        """
        Gets potentially rotatable bonds in molecule

        First, it retrieves all the atoms. Then, the bonds are filtered,
        discarding coordination (pseudo)bonds and sort them by atom serial.

        For each bond, try to retrieve it from the cache. If unavailable,
        request a bond rotation object to chimera.BondRot.

        In this step, we have to discard non rotatable atoms (as requested
        by the user), and make sure the involved atoms are of compatible type.
        Namely, one of them must be either C3, N3, C2 or N2, and both of them, 
        non-terminal (more than one neighbor).

        If the bond is valid, get the BondRot object. Chimera will complain
        if we already have requested that bond previously, or if the bond is in a
        cycle. Handle those exceptions silently, and get the next bond in that case.

        If no exceptions are raised, then store the rotation anchor in the BondRot
        object (that's the nearest atom in the bond to the molecular anchor),
        and store the BondRot object in the rotations cache.
        """
        atoms = self.parent.genes[self.target].compound.mol.atoms
        bonds = set(b for a in atoms for b in a.bonds if not a.element.isMetal)
        bonds = sorted(
            bonds, key=lambda b: min(y.serialNumber for y in b.atoms))

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
        """
        Probably unneded now
        """
        self.rotatable_bonds[:] = list(self.get_rotatable_bonds())

    @property
    def anchor(self):
        """
        Get the molecular anchor. Ie, the *root* of the rotations, the fixed
        atom of the molecule.

        Usually, this is the target atom in the Search gene, but if we can't find it,
        get the ``donor`` atom of the molecule.
        """
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
        return anchor
