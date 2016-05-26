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
This objective calculates the angle formed by three
given atoms (or the dihedral, if four atoms are given) and returns
the absolute difference of that angle and the target value.

"""

# Python
from __future__ import print_function
import math
import logging
# Chimera
import chimera
# GAUDI
from gaudi.objectives import ObjectiveProvider
from gaudi import parse

logger = logging.getLogger(__name__)


def enable(**kwargs):
    kwargs = Angle.validate(kwargs)
    return Angle(**kwargs)


class Angle(ObjectiveProvider):

    """
    Angle class

    Parameters
    ----------
    threshold : float
        Optimum angle
    tolerance :
        Allowed difference
    probes : list of str
        Atoms that make the angle, expressed as a series of 
        <molecule_name>/<serial_number> strings
    """

    validate = parse.Schema({
        parse.Required('probes'): parse.AssertList(parse.Named_spec("molecule", "atom")),
        parse.Required('threshold'): parse.Any(parse.Coerce(float), parse.In(['planar'])),
        'tolerance': parse.Coerce(float),
        }, extra=parse.ALLOW_EXTRA)

    def __init__(self, threshold=None, tolerance=-0.1, probes=None,
                 *args, **kwargs):
        ObjectiveProvider.__init__(self, **kwargs)
        self.threshold = threshold
        self.tolerance = tolerance
        self._probes = probes

    def probes(self, ind):
        for probe in self._probes:
            mol, serial = probe
            try:
                if isinstance(serial, int):
                    atom = next(a for a in ind.genes[mol].compound.mol.atoms
                                if serial == a.serialNumber)
                else:
                    atom = next(a for a in ind.genes[mol].compound.mol.atoms
                                if serial == a.name)
            except KeyError:
                logger.exception("Molecule not found")
                raise
            except StopIteration:
                logger.exception("No atoms matched for probe %s", probe)
                raise
            else:
                yield atom

    def evaluate(self, ind):
        atoms_coords = [a.molecule.openState.xform.apply(a.coord()) for a in self.probes(ind)]
        delta = 180.0
        try:
            angle = chimera.angle(*atoms_coords)
        except TypeError:  # four atoms, means dihedral
            angle = chimera.dihedral(*atoms_coords)
        except TypeError:  # threshold is str, calc abs sine
            if self.threshold == 'planar':
                delta = abs(math.sin(math.radians(angle)))
        else:
            delta = self.threshold - angle.real

        return delta

# TODO: Probes get lost if rotamers are applied!
