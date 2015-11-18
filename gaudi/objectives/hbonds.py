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
This objective is a wrapper around Chimera's `FindHBond`. 
It returns the number of hydrogen bonds that can be formed
between the target molecule and its environment.

.. todo::

    Evaluate the possible HBonds with some kind of function that
    gives a rough idea of the strength (energy) of each of them.

"""

# Python
import logging
# Chimera
import chimera
import FindHBond
# GAUDI
from gaudi.objectives import ObjectiveProvider
import gaudi.box

logger = logging.getLogger(__name__)


def enable(**kwargs):
    return Hbonds(**kwargs)


class Hbonds(ObjectiveProvider):

    """
    Hbonds class

    Parameters
    ----------
    probe : str
        Name of molecule being object of analysis
    radius : float
        Maximum distance from any point of probe that is searched
        for a possible interaction
    distance_tolerance : float, optional
        Allowed deviation from ideal distance to consider a valid H bond.
    angle_tolerance : float, optional
        Allowed deviation from ideal angle to consider a valid H bond.

    """

    def __init__(self, probe=None, radius=5.0, distance_tolerance=0.4, angle_tolerance=20.0,
                 *args, **kwargs):
        ObjectiveProvider.__init__(self, **kwargs)
        self._probe = probe
        self.distance_tolerance = distance_tolerance
        self.angle_tolerance = angle_tolerance
        self.radius = radius

    def molecules(self, ind):
        return tuple(m.compound.mol for m in ind.genes.values()
                     if m.__class__.__name__ == "Molecule")

    def probe(self, ind):
        return ind.genes[self._probe].compound.mol

    def evaluate(self, ind):
        """
        Find H bonds within self.radius angstroms from self.probe, and return
        only those that interact with probe. Ie, discard those hbonds in that search
        space whose none of their atoms involved are not part of self.probe.
        """
        molecules = self.molecules(ind)
        probe = self.probe(ind)
        test_atoms = self.surrounding_atoms(probe, molecules)
        hbonds = FindHBond.findHBonds(molecules, cacheDA=self._cache,
                                      donors=test_atoms, acceptors=test_atoms,
                                      distSlop=self.distance_tolerance,
                                      angleSlop=self.angle_tolerance)
        hbonds = FindHBond.base.filterHBondsBySel(
            hbonds, probe.atoms, 'any')

        return len(hbonds)

    def display(self, bonds):
        """
        Mock method to show a graphical depiction of the found H Bonds.
        """
        return gaudi.box.draw_interactions(bonds, name=self.name,
                                           startCol='00FFFF', endCol='00FFFF')

    ###
    def surrounding_atoms(self, probe, molecules):
        self.zone.clear()
        self.zone.add(probe.atoms)
        self.zone.merge(chimera.selection.REPLACE,
                        chimera.specifier.zone(
                            self.zone, 'atom',
                            None, self.radius,
                            molecules
                        )
                        )
        return self.zone.atoms()
