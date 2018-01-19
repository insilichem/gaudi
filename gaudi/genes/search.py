#!/usr/bin/env python
# -*- coding: utf-8 -*-

##############
# GaudiMM: Genetic Algorithms with Unrestricted
# Descriptors for Intuitive Molecular Modeling
#
# https://github.com/insilichem/gaudi
#
# Copyright 2017 Jaime Rodriguez-Guerra, Jean-Didier Marechal
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
##############

"""
This module provides spatial exploration of the environment.

It works by creating a sphere with radius `self.radius` and origin at
`self.origin`. The movement is achieved with three matrices that contain
a translation, a rotation, and a reference position.

It depends on :class:`gaudi.genes.molecule.Molecule`, since these are the ones
that will be moved around. Combined with the adequate objectives, this module
can be used to implement docking experiments.

"""

# Python
import random
import logging
from numpy import around as numpy_around
# Chimera
import chimera
from chimera import Xform as X
import Matrix as M
from FitMap.search import random_rotation
# GAUDI
from gaudi.genes import GeneProvider
from gaudi import parse


ZERO = chimera.Point(0.0, 0.0, 0.0)
IDENTITY = ((1.0, 0.0, 0.0, 0.0),
            (0.0, 1.0, 0.0, 0.0),
            (0.0, 0.0, 1.0, 0.0))
logger = logging.getLogger(__name__)


def enable(**kwargs):
    kwargs = Search.validate(kwargs)
    return Search(**kwargs)


class Search(GeneProvider):

    """
    Parameters
    ----------
    target : namedtuple
        The *anchor* atom of the molecule we want to move, with syntax
        ``<molecule_name>/<index>``. For example, if we want to move Ligand
        using atom with serial number = 1 as pivot, we would specify
        ``Ligand/1``. It's parsed to the actual chimera.Atom later on.
    center : 3-item list or tuple of float, optional
        Coordinates to the center of the desired search sphere
    radius : float
        Maximum distance from center that the molecule can move
    rotate : bool, bool
        If False, don't rotatate the molecule - only translation
    precision : int, bool
        Rounds the decimal part of the 3D search matrix to get a coarser
        model of space. Ie, less points can be accessed, the search is less
        exhaustive, more variability in less runs.

    Attributes
    ----------
    allele : 3-tuple of 4-tuple of floats
        A 4x3 matrix of float, as explained in Notes.
    origin : 3-tuple of float
        The initial position of the requested target molecule. If we don't take this
        into account, we can't move the molecule around was not originally in the
        center of the sphere.

    Notes
    -----
    **How matricial translation and rotation takes place**

    A single movement is summed up in a 4x3 matrix:

        (
        (R1, R2, R3, T1),
        (R4, R5, R6, T2),
        (R7, R8, R9, T3)
        )

    R-elements contain the rotation information, while T elements account for
    the translation movement.

    That matrix can be obtained from multipying three different matrices with
    this expression:

        multiply_matrices(translation, rotation, to_zero)

    To understand the operation, it must be read from the right:

        1. First, translate the molecule the origin of coordinates 0,0,0
        2. In that position, the rotation can take place.
        3. Then, translate to the final coordinates from zero. There's no need
           to get back to the original position.

    How do we get the needed matrices?

    - ``to_zero``. Record the original position (`origin`) of the molecule and
      multiply it by -1. Done with method `to_zero()`.

    - ``rotation``. Obtained directly from ``FitMap.search.random_rotation``

    - ``translation``. Check docstring of ``random_translation()`` in this module.

    """

    _validate = {
        parse.Required('target'): parse.Named_spec("molecule", "atom"),
        'center': parse.Any(parse.Coordinates, parse.Named_spec("molecule", "atom")),
        'radius': parse.Coerce(float),
        'rotate': parse.Boolean,
        'precision': parse.All(parse.Coerce(int), parse.Range(min=-3, max=6)),
        'interpolation': parse.All(parse.Coerce(float), parse.Range(min=0, max=1.0))
        }

    def __init__(self, target=None, center=None, radius=None, rotate=True,
                 precision=0, interpolation=0.5, **kwargs):
        GeneProvider.__init__(self, **kwargs)
        self.radius = radius
        self.rotate = rotate
        self.precision = precision
        self._center = center
        self.target = target
        self.interpolation = interpolation

    def __ready__(self):
        self.allele = self.random_transform()

    @property
    def center(self):
        if self._center:
            return parse_origin(self._center, self.parent)
        else:
            return self.origin

    @property
    def molecule(self):
        return self.parent.find_molecule(self.target.molecule).compound.mol

    @property
    def origin(self):
        return parse_origin(self.target, self.parent)

    @property
    def to_zero(self):
        """
        Return a translation matrix that takes the molecule from its
        original position to the origin of coordinates (0,0,0).

        Needed for rotations.
        """
        x, y, z = self.origin
        return ((1.0, 0.0, 0.0, -x),
                (0.0, 1.0, 0.0, -y),
                (0.0, 0.0, 1.0, -z))

    def express(self):
        """
        Multiply all the matrices, convert the result to a chimera.CoordFrame and
        set that as the xform for the target molecule. If precision is set, round them.
        """
        matrices = self.allele + (self.to_zero,)
        if self.precision > 0:
            self.molecule.openState.xform = M.chimera_xform(
                M.multiply_matrices(*numpy_around(matrices, self.precision).tolist()))
        else:
            self.molecule.openState.xform = M.chimera_xform(
                M.multiply_matrices(*matrices))

    def unexpress(self):
        """
        Reset xform to unity matrix.
        """
        self.molecule.openState.xform = X()

    def mate(self, mate):
        """
        Interpolate the matrices and assign them to each individual.
        Ind1 gets the rotated interpolation, while Ind2 gets the translation.
        """
        xf1 = M.chimera_xform(M.multiply_matrices(*self.allele))
        xf2 = M.chimera_xform(M.multiply_matrices(*mate.allele))
        interp = M.xform_matrix(M.interpolate_xforms(xf1, ZERO, xf2, 0.5))
        interp_rot = [x[:3] + (0,) for x in interp]
        interp_tl = [y[:3] + x[-1:] for x, y in zip(interp, M.identity_matrix())]
        self.allele, mate.allele = (self.allele[0], interp_rot), (interp_tl, mate.allele[1])

    def mutate(self, indpb):
        if random.random() < self.indpb:
            xf1 = M.chimera_xform(M.multiply_matrices(*self.allele))
            xf2 = M.chimera_xform(M.multiply_matrices(*self.random_transform()))
            interp = M.xform_matrix(M.interpolate_xforms(xf1, ZERO, xf2, self.interpolation))
            interp_rot = [x[:3] + (0,) for x in interp]
            interp_tl = [y[:3] + x[-1:] for x, y in zip(interp, M.identity_matrix())]
            self.allele = interp_tl, interp_rot

    #####
    def random_transform(self):
        """
        Wrapper function to provide translation and rotation in a single call
        """
        rotation = random_rotation() if self.rotate else IDENTITY
        translation = random_translation(self.center, self.radius)
        return translation, rotation


#############
# Some useful functions
def translate(molecule, anchor, target):
    if isinstance(anchor, chimera.Atom):
        anchor = anchor.coord()
    if isinstance(target, chimera.Atom):
        target = target.coord()

    t = X.translation(target - anchor)
    for a in molecule.atoms:
        a.setCoord(t.apply(a.coord()))


def rotate(molecule, at, alpha):
    if len(at) == 3:
        try:
            a1, a2, a3 = [a.coord() for a in at]
        except AttributeError:
            a1, a2, a3 = at
        axis_a = a1 - a2
        axis_b = a3 - a2
        delta = chimera.angle(a1, a2, a3) - alpha
        axis = chimera.cross(axis_a, axis_b)
        if axis.data() == (0.0, 0.0, 0.0):
            axis = chimera.cross(axis_a, axis_b + chimera.Vector(1, 0, 0))
            logger.warning("Had to choose arbitrary normal vector")
        pivot = a2
    elif len(at) == 4:
        try:
            a1, a2, a3, a4 = [a.coord() for a in at]
        except AttributeError:
            a1, a2, a3, a4 = at
        axis = a3 - a2
        delta = chimera.dihedral(a1, a2, a3, a4) - alpha
        pivot = a3
    else:
        raise ValueError(
            "Atom list must contain 3 (angle) or 4 (dihedral) atoms only")

    r = X.translation(pivot - ZERO)  # move to origin
    r.multiply(X.rotation(axis, - delta))  # rotate
    r.multiply(X.translation(ZERO - pivot))  # return to orig pos
    for a in molecule.atoms:
        a.setCoord(r.apply(a.coord()))


def rand_xform(origin, destination, r, rotate=True):
    x, y, z = origin
    to_zero = ((1.0, 0.0, 0.0, -x),
               (0.0, 1.0, 0.0, -y),
               (0.0, 0.0, 1.0, -z))
    rotation = random_rotation() if rotate else IDENTITY
    translation = random_translation(destination, r)
    return translation, rotation, to_zero


def random_translation(center, r):
    """
    Get a random point from the cube built with l=r and test if it's within
    the sphere. Most of the points will be, but not all of them, so get another
    one until that criteria is met.

    Parameters
    ----------
    center : 3-tuple of float
        Coordinates of the center of the search sphere
    r : float
        Radius of the search sphere

    Returns
    -------
    A translation matrix to a random point in the search sphere, with no rotation.
    """
    r2 = r*r
    a, b, c = center
    random_uniform = random.uniform
    while True:
        x, y, z = [random_uniform(-r, r) for m in center]
        if x*x + y*y + z*z <= r2:
            break
    return ((1.0, 0.0, 0.0, a + x),
            (0.0, 1.0, 0.0, b + y),
            (0.0, 0.0, 1.0, c + z))


def parse_origin(origin, individual=None):
    """
    The center of the sphere can be given as an Atom, or directly as
    a list of three floats (x,y,z). If it's an Atom, find it and return
    the xyz coords. If not, just turn the list into a tuple.

    Parameters
    ----------
    origin : 3-item list of coordinates, or chimera.Atom
    genes : gaudi.parse.Settings.genes
        List of gene-dicts to look for molecules that may contain
        the referred atom in `origin`

    Returns
    -------
    Tuple of float
        The x,y,z coordinates
    """
    if individual and isinstance(origin, tuple) and len(origin) == 2:
        mol, serial = origin
        return individual.find_molecule(mol).find_atom(serial).coord().data()
    elif isinstance(origin, list) and len(origin) == 3:
        return tuple(origin)
    else:
        raise ValueError('Origin {} cannot be parsed'.format(origin))
