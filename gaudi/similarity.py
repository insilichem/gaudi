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
This module contains the similarity functions that are used
to discard individuals that are not different enough.

"""
from __future__ import print_function, division
import logging

logger = logging.getLogger(__name__)


def rmsd(ind1, ind2, subjects, threshold, *args, **kwargs):
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
    for s in subjects:
        if s not in ind1.genes:
            raise ValueError('Molecule {} not found in individual'.format(s))

    molecules1 = [ind1._molecules[s] for s in subjects]
    molecules2 = [ind2._molecules[s] for s in subjects]

    # If ligands are not the same molecule, of course they aren't similar
    alleles1 = [g.allele for g in molecules1]
    alleles2 = [g.allele for g in molecules2]
    if alleles1 != alleles2:
        return False

    logger.debug("Comparing RMSD between #%s and #%s", id(ind1), id(ind2))
    rmsds = []
    for m1, m2 in zip(molecules1, molecules2):
        coords1 = m1._expressed_xformcoords
        coords2 = m2._expressed_xformcoords
        if coords1.shape[0] != coords2.shape[0]:
            return False
        rmsd_squared = ((coords1-coords2)**2).sum() / coords1.shape[0]
        rmsds.append(rmsd_squared)
    logger.debug("RMSD: " + str(rmsds))
    return all(rmsd < threshold*threshold for rmsd in rmsds)
