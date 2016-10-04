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
    """
    _validate = {
        parse.Required('probes'): parse.AssertList(parse.Named_spec("molecule", "atom")),
        parse.Required('target'): parse.Named_spec("molecule", "atom"),
        parse.Required('threshold'): parse.Any(parse.Coerce(float), parse.In(['covalent'])),
        'tolerance': parse.Coerce(float)
    }

    def __init__(self, threshold=None, tolerance=None, target=None, probes=None,
                 *args, **kwargs):
        ObjectiveProvider.__init__(self, **kwargs)
        self.threshold = threshold
        self.tolerance = tolerance
        self._probes = probes
        self._target = target

    def atoms(self, ind, *targets):
        for target in targets:
            mol, serial = target
            try:
                atom = next(a for a in ind.genes[mol].compound.mol.atoms
                            if serial == a.serialNumber)
            except StopIteration:
                logger.exception("No atoms matched for target %s", target)
                raise
            else:
                yield atom

    def evaluate(self, ind):
        """
        Measure the distance
        """
        distances = []
        target = next(self.atoms(ind, self._target))
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

    @staticmethod
    def _distance(atom1, atom2):
        return atom1.xformCoord().distance(atom2.xformCoord())
