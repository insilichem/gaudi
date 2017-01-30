#!/usr/bin/env python
# -*- coding: utf-8 -*-

##############
# GaudiMM: Genetic Algorithms with Unrestricted
# Descriptors for Intuitive Molecular Modeling
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
    molecules1 = [ind1.find_molecule(s) for s in subjects]
    molecules2 = [ind2.find_molecule(s) for s in subjects]

    logger.debug("Comparing RMSD between #%s and #%s", id(ind1), id(ind2))
    for m1, m2 in zip(molecules1, molecules2):
        coords1 = m1._expressed_coordinates
        coords2 = m2._expressed_coordinates
        if coords1.shape[0] != coords2.shape[0]:
            return False
        rmsd_squared = _rmsd_squared(coords1, coords2)
        if rmsd_squared > threshold*threshold:
            return False
    return True


def _rmsd_squared(coords1, coords2):
    diff = coords1 - coords2
    return (diff * diff).sum() / coords1.shape[0]
