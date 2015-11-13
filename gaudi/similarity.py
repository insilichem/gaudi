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
:mod:`gaudi.base` contains the core classes we use to build individuals
(potential solutions of the optimization process).
"""

import logging
from math import sqrt

logger = logging.getLogger(__name__)


def rmsd(ind1, ind2, subject, threshold):
    """
    Returns the RMSD between two individuals

    Parameters
    ----------
    ind1, ind2 : gaudi.base.Individual
    subject : str
        Name of gaudi.genes.molecule instance to measure
    threshold : float
        Maximum RMSD value to consider two individuals as similar.
        If rmsd > threshold, they are considered different.

    Returns
    -------
    bool
        True if rmsd is within threshold, False otherwise
    """
    logger.debug("Comparing RMSD between #%s and #%s",
                 id(ind1), id(ind2))
    ind1.express()
    compound1 = next(g for g in ind1.genes.values()
                     if g.__class__.__name__ == 'Molecule'
                     and g.name == subject).compound
    atoms1 = sorted(compound1.mol.atoms, key=lambda x: x.serialNumber)
    coords1 = [a.coord() for a in atoms1]
    xf1 = compound1.mol.openState.xform
    ind1.unexpress()

    ind2.express()
    compound2 = next(g for g in ind2.genes.values()
                     if g.__class__.__name__ == 'Molecule'
                     and g.name == subject).compound
    atoms2 = sorted(compound2.mol.atoms, key=lambda x: x.serialNumber)
    coords2 = [a.coord() for a in atoms2]
    xf2 = compound2.mol.openState.xform
    ind2.unexpress()

    sqdist = sum(xf1.apply(a).sqdistance(xf2.apply(b))
                 for a, b in zip(coords1, coords2))
    rmsd = sqrt(sqdist / ((len(coords1) + len(coords2)) / 2.0))
    logger.debug("RMSD: %f", rmsd)
    return rmsd < threshold
