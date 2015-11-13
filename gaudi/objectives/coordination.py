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
This objective performs rough estimations
of good orientations of ligating residues in a protein to
coordinate a given metal or small molecule. The geometry is approximated
by computing average distances from ligating atoms the metal centre (`self.probe`)
as well as the angles formed by the probe, the ligating atom and its immediate neighbor.
Good planarity is assured by a dihedral check.

"""

# Python
import math
import logging
import numpy
# Chimera
import chimera
import MetalGeom
# GAUDI
from gaudi.objectives import ObjectiveProvider
import gaudi.parse

logger = logging.getLogger(__name__)


def enable(**kwargs):
    return SimpleCoordination(**kwargs)


class SimpleCoordination(ObjectiveProvider):

    """
    SimpleCoordination class

    Parameters
    ----------
    probe : str
        The atom that acts as the metal center, expressed as 
        <molecule_name>/<atom serial>. This will be parsed later on.
    radius : float
        Distance from `probe` where ligating atoms must be found
    atom_types : list of str
        Types of atoms that are considered ligands to `probe`
    residues : list of str
        Type of residues that must coordinate to `probe`, expressed as
        <molecule_name>/<residue position>
    distance : float
        Target distance from `probe` to ligating atoms
    angle : float
        Target angle `probe`, ligand and neighbor should form ideally
    """

    def __init__(self, method='simple', probe=None, radius=None, atom_types=None, residues=None,
                 distance=None, angle=None, dihedral=None, min_atoms=1,
                 enforce_all_residues=False, only_one_ligand_per_residue=False, *args, **kwargs):
        ObjectiveProvider.__init__(self, **kwargs)
        self.method = method
        self._probe = probe
        self._residues = residues
        self.radius = radius
        self.atom_types = atom_types
        self.distance = distance
        self.angle = angle
        self.dihedral = dihedral
        self.min_atoms = min_atoms
        self.only_one_ligand_per_residue = only_one_ligand_per_residue
        self.enforce_all_residues = enforce_all_residues
        if self.method == 'metalgeom':
            self.evaluate = self.evaluate_MetalGeom
        else:
            self.evaluate = self.evaluate_simple

    def probe(self, ind):
        return self._getatom(ind, self._probe)

    def molecules(self, ind):
        return tuple(g.compound.mol for g in ind.genes.values()
                     if g.__class__.__name__ == "Molecule")

    def residues(self, ind):
        return set(self._getresidue(ind, *self._residues))

    def evaluate_simple(self, ind):
        """
        1. For every atom within the search radius that matches the criteria, get:
            - Coordinates for `self.probe`, atom, its immediate neighbor, and next neighbor
              not being terminal
            - Angle formed by probe, atom, and neighbor
            - Dihedral formed by probe, atom, neighbor, and next neighbor
        2. The score returned is the sum of averages of:
            - Absolute difference between `self.distance` and found atoms. That way,
              atoms that are within requested distance, get scores near zero.
            - Absolute difference of sines of `self.angle` and formed angles. Good matches
              tend to zero.
            - Absolute sine of dihedrals. If they are coplanar, this should be zero.
        3. As a result, lower scores are better.
        """
        try:
            atoms_by_distance = self._get_nearest_atoms()
        except (NotEnoughAtomsError, SomeResiduesMissingError):
            return -1000 * self.weight
        else:
            distances, angles, dihedrals = [], [], []
            for d, a2 in atoms_by_distance:
                # Get coords for needed atoms
                distances.append(d)
                if self.dihedral or self.angle:
                    c1 = self._get_xform_coord(self.probe)  # the probe
                    c2 = self._get_xform_coord(a2)  # the ligand is a2
                    a3 = a2.neighbors[0]  # bonded atom to ligand
                    c3 = self._get_xform_coord(a3)
                    a4 = next(a for a in a3.neighbors  # neighbor of a3 that is not ligand
                              if a is not a3 and len(a.neighbors) > 1)
                    c4 = self._get_xform_coord(a4)

                    if self.angle:
                        angle = chimera.angle(c1, c2, c3)
                        delta = abs(math.cos(math.radians(self.angle))
                                    - math.cos(math.radians(angle)))
                        angles.append(delta)

                    if self.dihedral:
                        dihedral = chimera.dihedral(c1, c2, c3, c4)
                        planarity = abs(math.sin(math.radians(dihedral)))
                        dihedrals.append(planarity)

            return sum(numpy.average(x) for x in (distances, angles, dihedrals) if x)

    def evaluate_MetalGeom(self, ind):
        """
        1. Get requested atoms sorted by distance
        2. If they meet the minimum quantity, return the rmsd for
           that geometry
        3. If that's not possible of they are not enough, return penalty
        """
        try:
            test_atoms = [a for d, a in self._get_nearest_atoms()]
        except (NotEnoughAtomsError, SomeResiduesMissingError):
            return -1000 * self.weight
        else:
            try:
                rmsd, center, vectors = MetalGeom.gui.geomDistEval(
                    self.geometry, self.target, test_atoms)
            except:  #
                logger.warning(
                    "Geometry not feasible in current conditions")
                rmsd = -1000 * self.weight
            return rmsd

    def _get_nearest_atoms(self):
        """
        1. Get atoms and residues found within `self.radius` angstroms from `self.probe`
        1.1. Found residues MUST include self.residues. Otherwise, apply penalty
        2. Sort atoms by absolute difference of `self.distance` and distance to `self.probe`.
           That way, nearest atoms are computed first.
        2.1. If found atoms do not include some of the requested types, apply penalty.
        """
        self._update_zone()
        atoms = self.zone.atoms()
        # (distance, ligand) tuple, sorted by distances
        atoms_by_distance = \
            [(abs(self.distance - self.probe.xformCoord().distance(a.xformCoord())),
              a) for a in atoms if a.name in self.atom_types and a.residue in self.residues]
        if len(atoms_by_distance) < self.min_atoms:
            logger.warning(
                "Could not find requested atoms from residues in probe environment")
            raise NotEnoughAtomsError
        found_residues = set(a.residue for d, a in atoms_by_distance)
        if self.enforce_all_residues and found_residues != self.residues:
            logger.warning(
                "Some atoms found, but some residues are missing")
            raise SomeResiduesMissingError

        atoms_by_distance.sort()
        return atoms_by_distance

    # TODO: Probes get lost if rotamers are applied!
    def _getatom(self, ind, probe):
        """
        Parse `Molecule/serialNumber` string and return a chimera.Atom located at
        `Molecule` with serial number `serialNumber`
        """
        mol, serial = gaudi.parse.parse_rawstring(probe)
        try:
            if isinstance(serial, int):
                atom = next(a for a in ind.genes[mol].compound.mol.atoms
                            if serial == a.serialNumber)
            else:
                atom = next(a for a in ind.genes[mol].compound.mol.atoms
                            if serial == a.name)
        except KeyError:
            logger.exception("Molecule %s not found", mol)
            raise
        except StopIteration:
            logger.exception("No atoms matched for probe %s", probe)
            raise
        else:
            return atom

    def _getresidue(self, ind, *residues):
        """
        Parse `Molecule/position` string and return a chimera.Residue located at
        `Molecule` with serial number `position`
        """
        for r in residues:
            mol, pos = gaudi.parse.parse_rawstring(r)
            try:
                res = next(r for r in ind.genes[mol].compound.mol.residues
                           if pos == r.id.position)

            except KeyError:
                logger.exception("Molecule %s not found", mol)
                raise
            except StopIteration:
                logger.exception("No residues matched for pos %s", pos)
                raise
            else:
                yield res

    def _update_zone(self, ind):
        """
        Clear existing selection and add atoms within `self.radius` from
        `self.probe`, as long as they belong to one of `self.molecules`
        """
        self.zone.clear()
        self.zone.add(self.probe(ind))
        self.zone.merge(chimera.selection.REPLACE,
                        chimera.specifier.zone(
                            self.zone, 'atom', None, self.radius, self.molecules)
                        )
        return self.zone

    @staticmethod
    def _get_xform_coord(a):
        """
        Return transformed coordinates from atom `a`
        """
        return a.molecule.openState.xform.apply(a.coord())


class NotEnoughAtomsError(Exception):
    pass


class SomeResiduesMissingError(Exception):
    pass
