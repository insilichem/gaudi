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
Document this!
"""

# Python
import math
import logging
import numpy
# Chimera
import chimera
# GAUDI
from gaudi.objectives import ObjectiveProvider
import gaudi.parse

logger = logging.getLogger(__name__)


def enable(**kwargs):
    return SimpleCoordination(**kwargs)


class SimpleCoordination(ObjectiveProvider):

    def __init__(self, probe=None, radius=None, atomtypes=None, residues=None,
                 distance=None, angle=None, dihedral=None,
                 *args, **kwargs):
        ObjectiveProvider.__init__(self, **kwargs)
        self.probe = self._getatom(probe)
        self.residues = set(self._getresidue(*residues))
        self.radius = radius
        self.atomtypes = atomtypes
        self.distance = distance
        self.angle = angle
        self.dihedral = dihedral
        self.molecules = tuple(g.compound.mol for g in self.parent.genes.values()
                               if g.__class__.__name__ == "Molecule")

    def evaluate(self):
        self._update_env()
        atoms, residues = self.env.atoms(), self.env.residues()
        if not self.residues.issubset(set(residues)):
            logger.warning(
                "Requested residues were not found in surroundings")
            return self.weight * -1000
        # Distance
        # (distance, ligand) tuple, sorted by distances
        dist_atoms = \
            sorted((abs(self.distance - self.probe.xformCoord().distance(a.xformCoord())),
                    a) for a in atoms if a.name in self.atomtypes)
        if not dist_atoms:
            logger.warning(
                "Could not find requested atoms in probe environment")
            return self.weight * -1000

        # Prevent several atoms per residue.
        # Only nearest is added, since distances is sorted
        residues = set()
        distances, angles, dihedrals = [], [], []
        for d, a2 in dist_atoms:
            if a2.residue not in residues:
                # Get coords for needed atoms
                c1 = self._get_xform_coord(self.probe)  # the probe
                c2 = self._get_xform_coord(a2)  # the ligand is a2
                a3 = a2.neighbors[0]  # bonded atom to ligand
                c3 = self._get_xform_coord(a3)
                a4 = next(a for a in a3.neighbors  # neighbor of a3 that is not ligand
                          if a is not a3 and len(a.neighbors) > 1)
                c4 = self._get_xform_coord(a4)

                angle = chimera.angle(c1, c2, c3)
                dihedral = chimera.dihedral(c1, c2, c3, c4)
                logger.debug("Distance %s, angle %s, dihedral %s",
                             d, angle, dihedral)
                delta = abs(math.cos(math.radians(self.angle))
                            - math.cos(math.radians(angle)))
                planarity = abs(math.sin(math.radians(dihedral)))
                logger.debug("Delta %s, planarity %s", delta, planarity)

                distances.append(abs(self.distance - d))
                angles.append(delta)
                dihedrals.append(planarity)
                residues.add(a2.residue)

        return sum(numpy.average(x) for x in (distances, angles, dihedrals))

    # TODO: Probes get lost if rotamers are applied!
    def _getatom(self, probe):
        mol, serial = gaudi.parse.parse_rawstring(probe)
        try:
            if isinstance(serial, int):
                atom = next(a for a in self.parent.genes[mol].compound.mol.atoms
                            if serial == a.serialNumber)
            else:
                atom = next(a for a in self.parent.genes[mol].compound.mol.atoms
                            if serial == a.name)
        except KeyError:
            logger.exception("Molecule %s not found", mol)
            raise
        except StopIteration:
            logger.exception("No atoms matched for probe %s", probe)
            raise
        else:
            return atom

    def _getresidue(self, *residues):
        for r in residues:
            mol, pos = gaudi.parse.parse_rawstring(r)
            try:
                res = next(r for r in self.parent.genes[mol].compound.mol.residues
                           if pos == r.id.position)

            except KeyError:
                logger.exception("Molecule %s not found", mol)
                raise
            except StopIteration:
                logger.exception("No residues matched for pos %s", pos)
                raise
            else:
                yield res

    def _update_env(self):
        self.env.clear()
        self.env.add(self.probe)
        self.env.merge(chimera.selection.REPLACE,
                       chimera.specifier.zone(
                           self.env, 'atom', None, self.radius, self.molecules)
                       )
        return self.env

    @staticmethod
    def _get_xform_coord(a):
        return a.molecule.openState.xform.apply(a.coord())
