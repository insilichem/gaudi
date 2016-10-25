#!/usr/bin/env python
# -*- coding: utf-8 -*-

##############
# GAUDIasm: Genetic Algorithms for Universal
# Design Inference and Atomic Scale Modeling
# Authors:  Jaime Rodriguez-Guerra Pedregal
#            <jaime.rodriguezguerra@uab.cat>
#           Jean-Didier Marechal
#            <jeandidier.marechal@uab.cat>
# Web: https://bitbucket.org/insilichem/gaudi
##############

"""
This objective calculates the distance between two
given atoms. It returns the absolute difference between the calculated
distance and the target value.

"""

# Python
import numpy
import logging
# Chimera
import chimera
# GAUDI
from gaudi.objectives import ObjectiveProvider
from gaudi import parse

logger = logging.getLogger(__name__)


def enable(**kwargs):
    kwargs = Distance.validate(kwargs)
    return Distance(**kwargs)


class Distance(ObjectiveProvider):

    """
    Distance class

    Parameters
    ----------
    threshold : float
        Optimum distance to meet
    tolerance : float
        Maximum deviation from threshold that is not penalized
    target : str
        The atom to measure the distance to, expressed as 
        <molecule name>/<atom serial>
    probes : list of str
        The atoms whose distance to `target` is being measured,
        expressed as <molecule name>/<atom serial>. If more than one
        is provided, the average of all of them is returned
    center_of_mass : bool

    """
    _validate = {
        parse.Required('probes'): parse.AssertList(parse.Named_spec("molecule", "atom")),
        parse.Required('target'): parse.Named_spec("molecule", "atom"),
        parse.Required('threshold'): parse.Any(parse.Coerce(float), parse.In(['covalent'])),
        'tolerance': parse.Coerce(float),
        'center_of_mass': parse.Coerce(float)
    }

    def __init__(self, threshold=None, tolerance=None, target=None, probes=None,
                 center_of_mass=False, *args, **kwargs):
        ObjectiveProvider.__init__(self, **kwargs)
        self.threshold = threshold
        self.tolerance = tolerance
        self.center_of_mass = center_of_mass
        self._probes = probes
        self._target = target
        if self.center_of_mass:
            self.evaluate = self.evaluate_center_of_mass
        else:
            self.evaluate = self.evaluate_distances

    def atoms(self, ind, *targets):
        for target in targets:
            mol, serial = target
            for atom in ind.find_molecule(mol).find_atoms(serial):
                yield atom

    def evaluate_distances(self, ind):
        """
        Measure the distance
        """
        distances = []
        target = ind.find_molecule(self._target.molecule).find_atom(self._target.atom)
        for a in self.atoms(ind, *self._probes):
            d = self._distance(a, target)
            if self.threshold == 'covalent':
                threshold = chimera.Element.bondLength(a.element, target.element)
            else:
                threshold = self.threshold
            d = d - threshold
            if self.tolerance is not None and d < self.tolerance:
                distances.append(-1000 * self.weight)
            else:
                distances.append(d)

        return numpy.mean(numpy.absolute(distances))

    def evaluate_center_of_mass(self, ind):
        target = ind.find_molecule(self._target.molecule).find_atom(self._target.atom)
        probes = list(self.atoms(ind, *self._probes))
        center_of_mass = self._center(*probes)
        
        return target.xformCoord().distance(chimera.Point(*center_of_mass))

    @staticmethod
    def _distance(atom1, atom2):
        return atom1.xformCoord().distance(atom2.xformCoord())

    @staticmethod
    def _center(*atoms):
        coords, masses = [], []
        for a in atoms:
            coords.append(a.xformCoord())
            masses.append(a.element.mass)

        return numpy.average(coords, axis=0, weights=masses)