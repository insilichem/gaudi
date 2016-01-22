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
This module contains the similarity functions that are used
to discard individuals that are not different enough.

"""
from __future__ import print_function, division
import logging

logger = logging.getLogger(__name__)


def rmsd(ind1, ind2, subjects, threshold):
    """
    Returns the RMSD between two individuals

    Parameters
    ----------
    ind1, ind2 : gaudi.base.Individual
    subjects : list of str
        Name of gaudi.genes.molecule instances to measure
    threshold : float
        Maximum RMSD value to consider two individuals as similar.
        If rmsd > threshold, they are considered different.

    Returns
    -------
    bool
        True if rmsd is within threshold, False otherwise

    """
    # If ligands are not the same molecule, of course they aren't similar
    molecules1 = [g.allele for g in ind1.genes if g.__class__.__name__ == "Molecule"]
    molecules2 = [g.allele for g in ind2.genes if g.__class__.__name__ == "Molecule"]
    if molecules1 != molecules2:
        return False

    logger.debug("Comparing RMSD between #%s and #%s", id(ind1), id(ind2))
    coords1 = _molecules_xform_coords_by_name(ind1, subjects)
    coords2 = _molecules_xform_coords_by_name(ind2, subjects)

    sqdist = sum(a.sqdistance(b) for a, b in zip(coords1, coords2))
    rmsd_squared = sqdist / ((len(coords1) + len(coords2)) / 2.0)
    logger.debug("RMSD: %f", rmsd_squared)
    return rmsd_squared < threshold*threshold


def _molecules_xform_coords_by_name(individual, subjects):
    individual.express()
    compounds = []
    for subject in subjects:
        for gene in individual.genes.values():
            if gene.__class__.__name__ == 'Molecule' and gene.name == subject:
                compounds.append(gene.compound)
    atoms = [a for compound in compounds
             for a in sorted(compound.mol.atoms, key=lambda x: x.serialNumber)]
    xform_coords = [a.xformCoord() for a in atoms]
    individual.unexpress()
    return xform_coords
