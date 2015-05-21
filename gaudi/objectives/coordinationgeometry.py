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
:mod:`gaudi.objectives.coordinationgeometry` is a wrapper around Chimera's
`MetalGeom`. It tries to fit the best coordination geometry in a given system,
returning the RMSD of the found conformation as compared to the ideal polyhedron.
"""

# Python
import logging
# Chimera
import chimera
import MetalGeom
# GAUDI
from gaudi.objectives import ObjectiveProvider
import gaudi.parse

logger = logging.getLogger(__name__)


def enable(**kwargs):
    return CoordinationGeometry(**kwargs)


class CoordinationGeometry(ObjectiveProvider):

    def __init__(self, target=None, geometry=None, min_atoms=None,
                 atom_types=None, radius=None, threshold=None,
                 *args, **kwargs):
        ObjectiveProvider.__init__(self, **kwargs)
        self.geometry = MetalGeom.geomData.geometries[geometry]
        self.min_atoms = min_atoms
        self.atom_types = atom_types
        self.radius = radius
        self.threshold = threshold
        self.molecules = tuple(m.compound.mol for m in self.parent.genes
                               if m.__class__.__name__ == "Molecule")
        mol, serial = gaudi.parse.parse_rawstring(target)
        try:
            if isinstance(serial, int):
                atom = next(a for a in self.parent.genes[mol].compound.mol.atoms
                            if serial == a.serialNumber)
            else:
                atom = next(a for a in self.parent.genes[mol].compound.mol.atoms
                            if serial == a.name)
        except KeyError:
            logger.exception("Molecule not found")
            raise
        except StopIteration:
            logger.exception("No atoms matched for target %s", target)
            raise
        else:
            self.target = atom

    def evaluate(self):
        test_atoms = [a for a in self._surrounding_atoms() if not a == self.target and
                      a.name in self.atom_types]

        if len(test_atoms) >= self.min_atoms:
            try:
                rmsd, center, vectors = MetalGeom.gui.geomDistEval(
                    self.geometry, self.target, test_atoms)
            except:  # geometry not feasible in current conditions
                rmsd = -1000 * self.weight
            return rmsd
        else:
            return -1000 * self.weight

    ###
    def _surrounding_atoms(self):
        self.env.clear()
        self.env.add(self.target)
        self.env.merge(chimera.selection.REPLACE,
                       chimera.specifier.zone(
                           self.env, 'atom', None, self.radius, self.molecules
                       )
                       )
        return self.env.atoms()
