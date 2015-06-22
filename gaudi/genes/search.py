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
:mod:`gaudi.genes.search` provides spatial exploration of the environment.

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
import gaudi.parse


ZERO = chimera.Point(0.0, 0.0, 0.0)
IDENTITY = ((1.0, 0.0, 0.0, 0.0),
            (0.0, 1.0, 0.0, 0.0),
            (0.0, 0.0, 1.0, 0.0))
logger = logging.getLogger(__name__)


def enable(**kwargs):
    return Search(**kwargs)


class Search(GeneProvider):

    def __init__(self, target=None, center=None, radius=None, rotate=True,
                 precision=None, **kwargs):
        GeneProvider.__init__(self, **kwargs)
        self.radius = radius
        self.rotate = rotate
        self.precision = precision
        self._center = center
        self._target = target
        self.target, self.anchor = target.split('/')

    def __ready__(self):
        self.allele = self.random_transform()

    @property
    def center(self):
        return parse_origin(self._center, self.parent.genes)

    @property
    def molecule(self):
        return self.parent.genes[self.target].compound.mol

    @property
    def origin(self):
        return parse_origin(self._target, self.parent.genes)

    @property
    def to_zero(self):
        origin = self.origin
        return ((1.0, 0.0, 0.0, -origin[0]),
                (0.0, 1.0, 0.0, -origin[1]),
                (0.0, 0.0, 1.0, -origin[2]))

    def express(self):
        """
        Multiply all the matrices, convert the result to a chimera.CoordFrame and
        set that as the xform for the target molecule. If precision is set, round them.
        """
        matrices = self.allele + (self.to_zero,)
        if self.precision is not None:
            self.molecule.openState.xform = M.chimera_xform(
                M.multiply_matrices(*numpy_around(matrices, self.precision).tolist()))
        else:
            self.molecule.openState.xform = M.chimera_xform(
                M.multiply_matrices(*matrices))

    def unexpress(self):
        """
        Reset xform to unity matrix
        """
        self.molecule.openState.xform = X()

    def mate(self, mate):
        """
        Interpolate the matrices and assign them to each individual.
        Ind1 gets the rotated interpolation, while Ind2 gets the translation.
        """
        xf1 = M.chimera_xform(M.multiply_matrices(*self.allele))
        xf2 = M.chimera_xform(M.multiply_matrices(*mate.allele))
        interp = M.xform_matrix(M.interpolate_xforms(xf1, ZERO,
                                                     xf2, 0.5))
        interp_rot = [x[:3] + (0,) for x in interp]
        interp_tl = [y[:3] + x[-1:]
                     for x, y in zip(interp, M.identity_matrix())]
        self.allele, mate.allele = (self.allele[0], interp_rot, self.allele[-1]), \
            (interp_tl, mate.allele[1], mate.allele[-1])

    def mutate(self, indpb):
        if random.random() < self.indpb:
            # Careful! Mutation generates a whole NEW position (similar to eta ~= 0)
            # TODO: We could use a eta param in mutation by interpolating original and
            # a new random xform with a given `frac` parameter
            self.allele = self.random_transform()

    #####
    def random_transform(self):
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
    to_zero = ((1.0, 0.0, 0.0, -origin[0]),
               (0.0, 1.0, 0.0, -origin[1]),
               (0.0, 0.0, 1.0, -origin[2]))
    rotation = random_rotation() if rotate else IDENTITY
    translation = random_translation(destination, r)
    return translation, rotation, to_zero


def random_translation(center, r):
    """
    Get a random point from the cube built with l=r and test if it's within
    the sphere. Most of the points will be, but not all of them, so get another
    one until that criteria is met.
    """
    inside = True
    while inside:
        x, y, z = [random.uniform(a - r, a + r) for a in center]
        if x * x + y * y + z * z > r:
            inside = False
    return ((1.0, 0.0, 0.0, x),
            (0.0, 1.0, 0.0, y),
            (0.0, 0.0, 1.0, z))


def parse_origin(origin, genes=None):
    """
    The center of the sphere can be given as an Atom, or directly as
    a list of three floats (x,y,z). If it's an Atom, find it and return
    the xyz coords. If not, just turn the list into a tuple
    """
    if isinstance(origin, str) and genes:
        mol, serial = gaudi.parse.parse_rawstring(origin)
        try:
            if isinstance(serial, int):
                atom = next(a for a in genes[mol].compound.mol.atoms
                            if serial == a.serialNumber)
            else:
                atom = next(a for a in genes[mol].compound.mol.atoms
                            if serial == a.name)
        except (KeyError, AttributeError, StopIteration):  # atom not found
            raise
        else:
            return atom.coord().data()
    elif isinstance(origin, list):
        return tuple(origin)
