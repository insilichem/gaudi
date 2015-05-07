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

# Chimera
import chimera
import FindHBond
# GAUDI
from gaudi.objectives import ObjectiveProvider
import gaudi.box


def enable(**kwargs):
    return Hbonds(**kwargs)


class Hbonds(ObjectiveProvider):

    def __init__(self, probe=None, distance_tolerance=0.4, angle_tolerance=20.0, radius=5.0,
                 *args, **kwargs):
        ObjectiveProvider.__init__(self, **kwargs)
        self.distance_tolerance = distance_tolerance
        self.angle_tolerance = angle_tolerance
        self.radius = radius
        self.molecules = tuple(m.compound.mol for m in self.parent.genes
                               if m.__class__.__name__ == "Molecule")
        try:
            self.probe = self.parent.genes[probe].compound.mol
        except KeyError:
            raise
        except AttributeError:
            raise

    def evaluate(self):
        test_atoms = self.surrounding_atoms()
        hbonds = FindHBond.findHBonds(self.molecules, cacheDA=self._cache,
                                      donors=test_atoms, acceptors=test_atoms,
                                      distSlop=self.distance_tolerance,
                                      angleSlop=self.angle_tolerance)
        hbonds = FindHBond.base.filterHBondsBySel(
            hbonds, self.probe.atoms, 'any')

        return len(hbonds)

    def display(self, bonds):
        return gaudi.box.draw_interactions(bonds, name=self.name,
                                           startCol='00FFFF', endCol='00FFFF')

    ###
    def surrounding_atoms(self):
        self.env.clear()
        self.env.add(self.probe.atoms)
        self.env.merge(chimera.selection.REPLACE,
                       chimera.specifier.zone(
                           self.env, 'atom', None, self.radius, self.molecules
                       )
                       )
        return self.env.atoms()
