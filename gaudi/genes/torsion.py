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
This module helps explore small molecules flexibility.

It does so by performing bond rotations in the selected `gaudi.genes.molecule.Molecule`
objects, if they exhibit free bond rotations.

"""

# Python
import random
import logging
# Chimera
import chimera
# External dependencies
from deap.tools import cxSimulatedBinaryBounded, mutPolynomialBounded
# GAUDI
from gaudi.genes import GeneProvider
from gaudi import box, parse

logger = logging.getLogger(__name__)


def enable(**kwargs):
    kwargs = Torsion.validate(kwargs)
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
    anchor : str
        Molecule/atom_serial_number of reference atom for torsions
    rotatable_atom_types : list of str
        Which type of atom types (as in chimera.Atom.idatmType) should rotate.
        Defaults to ('C3', 'N3', 'C2', 'N2', 'P').
    rotatable_atom_names : list of str
        Which type of atom names (as in chimera.Atom.name) should rotate.
        Defaults to ().

    Notes
    -----

    .. todo ::

        `max_bonds` should be automatically computed, based on ligand
        expected composition (careful with block-built ligands...)

    """

    _validate = {
        parse.Required('target'): parse.Molecule_name,
        'flexibility': parse.Degrees,
        'max_bonds': parse.All(parse.Coerce(int), parse.Range(min=0)),
        'anchor': parse.Named_spec("molecule", "atom"),
        'rotatable_atom_types': [basestring],
        'rotatable_atom_names': [basestring],
        'rotatable_elements': [basestring],
        'non_rotatable_bonds': [parse.All([parse.Named_spec("molecule", "atom")],
                                          parse.Length(min=2, max=2))],
        'precision': parse.All(parse.Coerce(int), parse.Range(min=0, max=3))
        }

    BONDS_ROTS = {}

    def __init__(self, target=None, flexibility=360.0, max_bonds=None, anchor=None,
                 rotatable_atom_types=('C3', 'N3', 'C2', 'N2', 'P'), 
                 rotatable_atom_names=(), rotatable_elements=(),
                 non_rotatable_bonds=(), precision=1, **kwargs):
        GeneProvider.__init__(self, **kwargs)
        self._kwargs = kwargs
        self.target = target
        self.flexibility = 360.0 if flexibility > 360 else flexibility
        self.max_bonds = max_bonds
        self.rotatable_atom_types = rotatable_atom_types
        self.rotatable_atom_names = rotatable_atom_names
        self.rotatable_elements = rotatable_elements
        self.non_rotatable_bonds = non_rotatable_bonds
        self.precision = precision
        self._anchor = anchor
        self.allele = [self.random_angle() for i in xrange(50)]

    def __expression_hooks__(self):
        if self.max_bonds is None:
            self.max_bonds = len(list(self.rotatable_bonds))
        self.allele = [self.random_angle() for i in xrange(self.max_bonds)]

    def express(self):
        """
        Apply rotations to rotatable bonds
        """
        for alpha, br in zip(self.allele, self.rotatable_bonds):
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
        self_allele, mate_allele = cxSimulatedBinaryBounded(
            self.allele, mate.allele, eta=self.cxeta,
            low=-0.5 * self.flexibility, up=0.5 * self.flexibility)
        self.allele[:] = [round(n, self.precision) for n in self_allele]
        mate.allele[:] = [round(n, self.precision) for n in mate_allele]

    def mutate(self, indpb):
        allele, = mutPolynomialBounded(self.allele,
                                       indpb=self.indpb, eta=self.mteta,
                                       low=-0.5 * self.flexibility,
                                       up=0.5 * self.flexibility)
        self.allele[:] = [round(n, self.precision) for n in allele]

    def clear_cache(self):
        GeneProvider.clear_cache()
        self.BONDS_ROTS.clear()

    #####
    @property
    def molecule(self):
        return self.parent.find_molecule(self.target).compound.mol

    def random_angle(self):
        """
        Returns a random angle within flexibility limits
        """
        return round(random.uniform(-0.5 * self.flexibility, 0.5 * self.flexibility),
                     self.precision)

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
        try:
            return self.molecule._rotatable_bonds
        except AttributeError:
            self.molecule._rotatable_bonds = list(self._compute_rotatable_bonds())
            return self.molecule._rotatable_bonds

    def _compute_rotatable_bonds(self):
        bonds = sorted(self.molecule.bonds, 
                       key=lambda b: min(y.serialNumber for y in b.atoms))

        non_rotatable_bonds = []
        for atom_a, atom_b in self.non_rotatable_bonds:
            a = self.parent.find_molecule(atom_a.molecule).find_atom(atom_a.atom)
            b = self.parent.find_molecule(atom_b.molecule).find_atom(atom_b.atom)
            bond = a.findBond(b)
            if bond:
                non_rotatable_bonds.append(bond)
            else:
                logger.warning('Atoms {} and {} are not bonded!'.format(a, b))

        def conditions(*atoms):
            for a in atoms:
                # Must be satisfied by at least one atom
                if a.numBonds == 1 or a.element.isMetal:
                    return False
            for a in atoms:
                if a.name == 'DUM' or \
                   a.idatmType in self.rotatable_atom_types or \
                   a.name in self.rotatable_atom_names or \
                   a.element.name in self.rotatable_elements:
                    return True

        for b in bonds:
            if b in non_rotatable_bonds:
                continue
            if conditions(*b.atoms):
                try:
                    br = chimera.BondRot(b)
                except (chimera.error, ValueError) as v:
                    if "cycle" in str(v) or "already used" in str(v):
                        continue  # discard bonds in cycles and used!
                    raise
                else:
                    br.rotanchor = box.find_nearest(self.anchor, b.atoms)
                    yield br

    @property
    def anchor(self):
        """
        Get the molecular anchor. Ie, the *root* of the rotations, the fixed
        atom of the molecule.

        Usually, this is the target atom in the Search gene, but if we can't find it,
        get the ``donor`` atom of the molecule.
        """
        try:
            return self.molecule._rotation_anchor
        except AttributeError:
            pass
        if self._anchor is not None:
            mol, atom = self._anchor
            try:
                molecule_gene = self.parent.find_molecule(mol)
                anchor = molecule_gene.find_atom(atom)
            except StopIteration:
                pass
            else:
                self.molecule._rotation_anchor = anchor
                return anchor
                
        target_gene = self.parent.find_molecule(self.target)
        try:
            search_gene = next(g for g in self.parent.genes.values()
                               if g.__class__.__name__ == 'Search'
                               and g.target == self.target)
        except StopIteration:
            anchor = target_gene.compound.donor
        else:
            try:
                anchor = target_gene.find_atom(search_gene.anchor)
            except (StopIteration, AttributeError):
                anchor = target_gene.compound.donor
        self.molecule._rotation_anchor = anchor
        return anchor
