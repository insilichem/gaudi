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
This module provides spatial exploration of the environment
through a MD trajectory file.
"""

# Python
import random
import logging
# 3rd party
import chimera
import mdtraj
# GAUDI
from gaudi.genes import GeneProvider
from gaudi import parse
from gaudi.objectives.energy import Energy

logger = logging.getLogger(__name__)


def enable(**kwargs):
    kwargs = Trajectory.validate(kwargs)
    return Trajectory(**kwargs)


class Trajectory(GeneProvider):

    """
    Parameters
    ----------
    target : str
        The Molecule that contains the topology of the trajectory.
    path : str
        Path to a MD trajectory file, as supported by mdtraj.
    max_frame : int
        Last frame of the trajectory that can be loaded.
    stride : int, optional
        Only load one in every `stride` frames
    preload : bool, optional
        Load the full trajectory in memory to accelerate expression.
        Not recommended for large files!

    Attributes
    ----------
    allele : The index of a frame in the MD trajectory.
    _traj : Alias to the frames cache
    """

    _validate = {
        parse.Required('target'): parse.Molecule_name,
        parse.Required('path'): parse.ExpandUserPathExists,
        parse.Required('max_frame'): parse.All(parse.Coerce(int), parse.Range(min=1)),
        'stride': parse.All(parse.Coerce(int), parse.Range(min=1)),
        'preload': bool,
        }

    def __init__(self, target=None, path=None, max_frame=None, stride=1, preload=False,
                 **kwargs):
        GeneProvider.__init__(self, **kwargs)
        self.target = target
        self.path = path
        self.max_frame = max_frame
        self.stride = stride
        self.preload = False
        try:
            self._traj = self._cache[self.name]
        except KeyError:
            self._traj = self._cache[self.name] = {}

    def __ready__(self):
        self.allele = self.random_frame_number()
        self._original_xyz = self.molecule.xyz(transformed=False)

    def __expression_hooks__(self):
        if self.preload and self.path not in self._traj:
            self._traj[self.path] = mdtraj.load(self.path, top=self.topology)

    @property
    def molecule(self):
        """
        The target Molecule gene
        """
        return self.parent.find_molecule(self.target)

    @property
    def topology(self):
        """
        Returns the equivalent mdtraj Topology object
        of the currently expressed Chimera molecule
        """
        mol = self.molecule.compound.mol
        try:
            return mol._mdtraj_topology
        except AttributeError:
            openmm_top = Energy.chimera_molecule_to_openmm_topology(mol)
            mdtraj_top = mdtraj.Topology.from_openmm(openmm_top)
            mol._mdtraj_topology = mdtraj_top
            return mdtraj_top

    def express(self):
        """
        Load the frame requested by the current allele into
        a new CoordSet object (always at index 1) and set 
        that as the active one.
        """
        traj = self.load_frame(self.allele)
        coords = traj.xyz[0] * 10
        for a, xyz in zip(self.molecule.compound.mol.atoms, coords):
            a.setCoord(chimera.Point(*xyz))

    def unexpress(self):
        """
        Set the original coordinates (stored at mol.coordSets[0]) 
        as the active ones.
        """
        for a, xyz in zip(self.molecule.compound.mol.atoms, self._original_xyz):
            a.setCoord(chimera.Point(*xyz))

    def mate(self, mate):
        """
        Simply exchange alleles. Can't try to interpolate 
        an intermediate structure because the result wouldn't
        probably belong to the original trajectory!
        """
        self.allele, mate.allele = mate.allele, self.allele

    def mutate(self, indpb):
        if random.random() < indpb:
            self.allele = self.random_frame_number()

    def random_frame_number(self):
        return random.choice(range(0, self.max_frame, self.stride))

    def load_frame(self, n):
        if self.preload:
            return self._traj[self.path][n]
        return mdtraj.load_frame(self.path, self.allele, top=self.topology)
