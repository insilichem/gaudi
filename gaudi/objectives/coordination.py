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
This objective performs rough estimations
of good orientations of ligating residues in a protein to
coordinate a given metal or small molecule. The geometry is approximated
by computing average distances from ligating atoms the metal centre (`self.probe`)
as well as the angles formed by the probe, the ligating atom and its immediate neighbor.
Good planarity is assured by a dihedral check.

"""

# Python
from __future__ import print_function, division
import math
import logging
import numpy
from math import acos, degrees
# Chimera
import chimera
from chimera import cross, Xform, Plane, angle, Point, Vector
from chimera.match import matchPositions
from MetalGeom.geomData import geometries as MG_geometries
# GAUDI
from gaudi.objectives import ObjectiveProvider
from gaudi import parse

logger = logging.getLogger(__name__)


def enable(**kwargs):
    kwargs = SimpleCoordination.validate(kwargs)
    return SimpleCoordination(**kwargs)


class SimpleCoordination(ObjectiveProvider):

    """
    SimpleCoordination class

    Parameters
    ----------
    probe : tuple
        The atom that acts as the metal center, expressed as
        <molecule_name>/<atom serial>. This will be parsed later on.
    radius : float
        Distance from `probe` where ligating atoms must be found
    atom_types : list of str
        Types of atoms that are considered ligands to `probe`
    residues : list of str
        Type of residues that must coordinate to `probe`, expressed as
        <molecule_name>/<residue position>
    distance : float
        Target distance from `probe` to ligating atoms
    angle : float
        Target angle `probe`, ligand and neighbor should form ideally
    """
    validate = parse.Schema({
        parse.Required('probe'): parse.Named_spec("molecule", "atom"),
        'radius': parse.Coerce(float),
        'atom_types': [str],
        'residues': [parse.Named_spec("molecule", "residue")],
        'distance': parse.All(parse.Coerce(float), parse.Range(min=0)),
        'angle': parse.Coerce(float),
        'min_atoms': parse.All(parse.Coerce(int), parse.Range(min=0)),
        'geometry': parse.In(MG_geometries.keys()),
        'enforce_all_residues': parse.Coerce(bool),
        'only_one_ligand_per_residue': parse.Coerce(bool),
        'method': parse.In(['simple', 'metalgeom', 'metalgeom_directional'])
        }, extra=parse.ALLOW_EXTRA)
    
    def __init__(self, method='simple', probe=None, radius=None, atom_types=(), residues=(),
                 distance=0, angle=None, dihedral=None, min_atoms=1, geometry='tetrahedral',
                 enforce_all_residues=False, only_one_ligand_per_residue=False, *args, **kwargs):
        ObjectiveProvider.__init__(self, **kwargs)
        self.method = method
        self._probe = probe
        self._residues = residues
        self.radius = radius
        self.atom_types = atom_types
        self.distance = distance
        self.angle = angle
        self.dihedral = dihedral
        self.min_atoms = min_atoms
        self.only_one_ligand_per_residue = only_one_ligand_per_residue
        self.enforce_all_residues = enforce_all_residues
        if self.method == 'metalgeom':
            self.evaluate = self.evaluate_MetalGeom
            self.geometry = MG_geometries[geometry]
        elif self.method == 'metalgeom_directional':
            self.evaluate = self.evaluate_MetalGeom_directional
            self.geometry = MG_geometries[geometry]
        else:
            self.evaluate = self.evaluate_simple

    def probe(self, ind):
        return self._getatom(ind, self._probe)

    def molecules(self, ind):
        return tuple(g.compound.mol for g in ind.genes.values()
                     if g.__class__.__name__ == "Molecule")

    def residues(self, ind):
        return set(self._getresidue(ind, *self._residues))

    def evaluate_simple(self, ind):
        """
        1. For every atom within the search radius that matches the criteria, get:
            - Coordinates for `self.probe`, atom, its immediate neighbor, and next neighbor
              not being terminal
            - Angle formed by probe, atom, and neighbor
            - Dihedral formed by probe, atom, neighbor, and next neighbor
        2. The score returned is the sum of averages of:
            - Absolute difference between `self.distance` and found atoms. That way,
              atoms that are within requested distance, get scores near zero.
            - Absolute difference of sines of `self.angle` and formed angles. Good matches
              tend to zero.
            - Absolute sine of dihedrals. If they are coplanar, this should be zero.
        3. As a result, lower scores are better.
        """
        try:
            atoms_by_distance = self._nearest_atoms(ind)
        except (NotEnoughAtomsError, SomeResiduesMissingError):
            return -1000 * self.weight
        else:
            distances, angles, dihedrals = [], [], []
            for d, a2 in atoms_by_distance:
                # Get coords for needed atoms
                distances.append(d)
                if self.dihedral or self.angle:
                    c1 = self._get_xform_coord(self.probe)  # the probe
                    c2 = self._get_xform_coord(a2)  # the ligand is a2
                    a3 = a2.neighbors[0]  # bonded atom to ligand
                    c3 = self._get_xform_coord(a3)
                    a4 = next(a for a in a3.neighbors  # neighbor of a3 that is not ligand
                              if a is not a3 and len(a.neighbors) > 1)
                    c4 = self._get_xform_coord(a4)

                    if self.angle:
                        angle = chimera.angle(c1, c2, c3)
                        delta = abs(math.cos(math.radians(self.angle))
                                    - math.cos(math.radians(angle)))
                        angles.append(delta)

                    if self.dihedral:
                        dihedral = chimera.dihedral(c1, c2, c3, c4)
                        planarity = abs(math.sin(math.radians(dihedral)))
                        dihedrals.append(planarity)

            return sum(numpy.average(x) for x in (distances, angles, dihedrals) if x)


    def evaluate_MetalGeom(self, ind):
        """
        1. Get requested atoms sorted by distance
        2. If they meet the minimum quantity, return the rmsd for
           that geometry
        3. If that's not possible of they are not enough, return penalty
        """
        try:
            test_atoms = [a for d, a in self._nearest_atoms(ind)]
        except (NotEnoughAtomsError, SomeResiduesMissingError):
            logger.warning("Not enough atoms or some residues missing")
            return -1000 * self.weight
        else:
            try:
                max_ligands = len(self.geometry.normVecs)
                rmsd = geomDistEval_patched(
                    self.geometry, self.probe(ind).xformCoord(),
                    [a.xformCoord() for a in test_atoms[:max_ligands]],
                    min_ligands=self.min_atoms)
            except Exception as e:
                logger.exception(e)  #
                logger.warning("Geometry not feasible in current conditions")
                rmsd = -1000 * self.weight
            return rmsd

    def evaluate_MetalGeom_directional(self, ind):
        """
        1. Get requested atoms sorted by distance
        2. If they meet the minimum quantity, return the rmsd for
           that geometry
        3. If that's not possible of they are not enough, return penalty
        """
        try:
            test_atoms = [a for d, a in self._nearest_atoms(ind)]
        except (NotEnoughAtomsError, SomeResiduesMissingError):
            logger.warning("Not enough atoms or some residues missing")
            return -1000 * self.weight
        
        max_ligands = len(self.geometry.normVecs)
        ligands = test_atoms[:max_ligands]
        ligand_coords = [a.xformCoord() for a in ligands]
        metal = self.probe(ind)
        metal_coord = metal.xformCoord()

        # rmsd
        try:
            rmsd = geomDistEval_patched(self.geometry, metal_coord, 
                                        ligand_coords, min_ligands=self.min_atoms)
        except Exception as e:
            logger.exception(e)  #
            logger.warning("Geometry not feasible in current conditions")
            return -1000 * self.weight
        
        # directionality
        ligand_objects = [LigandHelper(ligand, ligands, metal) for ligand in ligands]
        directionality = directional_evaluation(ligand_objects)

        return rmsd + directionality


    def _nearest_atoms(self, ind):
        """
        1. Get atoms and residues found within `self.radius` angstroms from `self.probe`
        1.1. Found residues MUST include self.residues. Otherwise, apply penalty
        2. Sort atoms by absolute difference of `self.distance` and distance to `self.probe`.
           That way, nearest atoms are computed first.
        2.1. If found atoms do not include some of the requested types, apply penalty.
        """

        self._update_zone(ind)
        atoms = self.zone.atoms()
        metal = self.probe(ind)
        residues = self.residues(ind)

        # (distance, ligand) tuple, sorted by distances
        def abs_distance(a): return abs(self.distance - metal.xformCoord().distance(a.xformCoord()))
        if self.atom_types:
            def atom_in_types(a): return a.name in self.atom_types 
        else:
            def atom_in_types(a): return True

        atoms_by_distance = [(abs_distance(a), a) for a in atoms 
                             if atom_in_types(a) and a.residue in residues]
        if len(atoms_by_distance) < self.min_atoms:
            logger.warning("Could not find requested atoms from residues in probe environment")
            raise NotEnoughAtomsError

        found_residues = set(a.residue for d, a in atoms_by_distance)
        if self.enforce_all_residues and found_residues != residues:
            logger.warning("Some atoms found, but some residues are missing")
            raise SomeResiduesMissingError

        atoms_by_distance.sort()
        return atoms_by_distance

    # TODO: Probes get lost if rotamers are applied!
    def _getatom(self, ind, probe):
        """
        Parse `Molecule/serialNumber` string and return a chimera.Atom located at
        `Molecule` with serial number `serialNumber`
        """
        mol, serial = probe
        try:
            if isinstance(serial, int):
                atom = next(a for a in ind.genes[mol].compound.mol.atoms
                            if serial == a.serialNumber)
            else:
                atom = next(a for a in ind.genes[mol].compound.mol.atoms
                            if serial == a.name)
        except KeyError:
            logger.exception("Molecule %s not found", mol)
            raise
        except StopIteration:
            logger.exception("No atoms matched for probe %s", probe)
            raise
        else:
            return atom

    def _getresidue(self, ind, *residues):
        """
        Parse `Molecule/position` string and return a chimera.Residue located at
        `Molecule` with serial number `position`
        """
        for mol, pos in residues:
            try:
                if pos == '*':
                    yield next(r for r in ind.genes[mol].compound.mol.residues)
                else:
                    yield next(r for r in ind.genes[mol].compound.mol.residues
                               if pos == r.id.position)
            except KeyError:
                logger.exception("Molecule %s not found", mol)
                raise
            except StopIteration:
                logger.exception("No residues matched for pos %s", pos)
                raise

    def _update_zone(self, ind):
        """
        Clear existing selection and add atoms within `self.radius` from
        `self.probe`, as long as they belong to one of `self.molecules`
        """
        self.zone.clear()
        self.zone.add(self.probe(ind))
        self.zone.merge(chimera.selection.REPLACE,
                        chimera.specifier.zone(
                            self.zone, 'atom', None, self.radius, self.molecules(ind)))
        return self.zone

    @staticmethod
    def _get_xform_coord(a):
        """
        Return transformed coordinates from atom `a`
        """
        return a.molecule.openState.xform.apply(a.coord())


class NotEnoughAtomsError(Exception):
    pass


class SomeResiduesMissingError(Exception):
    pass


def geomDistEval_patched(geom, metal_coord, ligands_coords, min_ligands=2):
    """
    Return RMSD of given geometry with best polyhedron built with metal_coord
    as center and ligands_coords as possible vertices.

    Parameters
    ----------
    geom : MetalGeom.Geometry.Geometry
        Geometry of desired coordination. Get it from MetalGeom.geomData.geometries.
    metal_coord : chimera.Point
        Coordinates of metal ion
    ligands_coords : list of chimera.Point
        Coordinates of each of the potential ligand atoms
    min_ligands : int, optional
        Number of ligands that must be coordinating
    Returns
    -------
    rmsd : float
        RMSD value of ideal polyhedron alignment with given coordinates
    """

    if len(ligands_coords) == 1:
        return (0.0, None, None)
    mc = metal_coord
    lcoords = ligands_coords
    lvecs = [lc - mc for lc in lcoords]

    for lv in lvecs:
        lv.normalize()
    # (1), (2): pick primary, secondary
    npv, nsv = lvecs[:2]
    nsvpt = Point(*nsv.data())
    origin = Point()
    # (3) for each vector in the geometry
    best_rmsd = [1000]
    for i1, nv1 in enumerate(geom.normVecs):
        # (3a) orient towards primary
        rotAxis = cross(nv1, npv)
        cos = nv1 * npv
        ang = degrees(acos(cos))
        if rotAxis.length > 0:
            xf1 = Xform.rotation(rotAxis, ang)
        elif cos > 0:
            xf1 = Xform.identity()
        else:
            xf1 = Xform.rotation(1.0, 0.0, 0.0, 180.0)
        plane1 = Plane(origin, cross(npv, nsv))
        # (3b) for each other vector in the geometry
        for i2, nv2 in enumerate(geom.normVecs):
            if i1 == i2:
                continue
            # (3b1) orient towards secondary
            xfnv2 = xf1.apply(nv2)
            # Vectors have no rounding tolerance for equality tests
            # if npv == xfnv2 or npv == -xfnv2:
            ang = angle(npv, xfnv2)
            if ang < 0.00001 or ang > 179.99999:
                xf2 = Xform.identity()
            else:
                plane2 = Plane(origin, cross(npv, xfnv2))
                ang = angle(plane1.normal, plane2.normal)
                if plane1.distance(Point(*xfnv2.data())) > 0.0:
                    angles = [0.0 - ang, 180.0 - ang]
                else:
                    angles = [ang, ang - 180.0]
                bestD = None
                for ang in angles:
                    xf = Xform.rotation(npv, ang)
                    d = nsvpt.sqdistance(Point(*xf.apply(xfnv2)))
                    if bestD == None or d < bestD:
                        bestD = d
                        xf2 = xf
            # (3b2) find correspondences between remaining
            #   ligands and other vectors, closest first
            xf = Xform.identity()
            xf.multiply(xf2)
            xf.multiply(xf1)
            correspondences = {npv: nv1, nsv: nv2}
            used = set([i1, i2])
            xfnvs = [xf.apply(nv) for nv in geom.normVecs]
            diffs = []
            for i, xfnv in enumerate(xfnvs):
                if i in used:
                    continue
                for lv in lvecs[2:]:
                    diffs.append((angle(xfnv, lv), lv, i))
            diffs.sort()
            while len(used) < len(ligands_coords):
                ang, lv, i = diffs.pop(0)
                if i in used or lv in correspondences:
                    continue
                correspondences[lv] = geom.normVecs[i]
                used.add(i)
            if len(used) < min_ligands:
                continue
            # (3b3) get RMSD based on these pairings (incl. metal)
            realPts = [mc] + lcoords
            idealVecs = []
            for lv, lc in zip(lvecs, lcoords):
                nv = correspondences[lv]
                idealVecs.append(xf.apply(nv) * (lc - mc).length)
            idealPts = [origin] + [origin + iv for iv in idealVecs]
            rmsdXf, rmsd = matchPositions(realPts, idealPts)
            best_rmsd.append(rmsd)
    return min(best_rmsd)

def directional_evaluation(ligands_list):
    """
    Given a ligand get its bonding direction and compare it with the real one.

    Parameters
    ----------
    ligands_list : list of Ligand objects
        It contains the ligand candidates information.

    Returns
    -------
    sin_sum/len(ligands_list): float
        This value is taken as a scorer for the directional evaluation of the
        metal coordination.
    """
    sin_sum = 0
    for ligand in ligands_list:
        direction = bond_direction(ligand)
        if all(not v for v in direction):
            continue
        if ligand.bidentate:
            metal_bond = ligand.bond_dir_origin - ligand.metal_center.coord()
            linearity = abs(math.sin(math.radians(angle(direction, metal_bond))))
        else:
            linearity = abs(math.sin(math.radians(angle(direction, ligand.vector))))
        sin_sum += linearity
    return sin_sum/(len(ligands_list))


def bond_direction(ligand):
    """
    Obtains a preferential bonding direction. The directional evaluation is
    calculated according to this direction.

    Parameters
    ----------
    ligand : Ligand object
        It is the ligand candidate whose bonding direction, if any, will be
        obtained.

    Returns
    -------
    vec : chimera.Vector
        It is the vector whose direction matches with the ligand's bonding
        direction.
    """
    p0 = ligand.coord
    if ligand.terminal:
        # For aspartic & glutamic acids (bidentates)
        if ligand.bidentate:
            p1 = ligand.bond_dir_origin
            p2 = ligand.bidentate_mate.coord()
            vec1 = Vector(*(p0-p1).data())
            vec2 = Vector(*(p2-p1).data())
            rot_axis = cross(vec1, vec2)
            ang = angle(vec1, vec2)/2
            rotation = Xform.rotation(rot_axis, ang)
            vec = rotation.apply(vec1)
            return vec

    # For histidine & tryptophan
    if len(ligand.neighbors) == 2:
        p1 = ligand.neighbors[0].coord()
        p2 = ligand.neighbors[1].coord()
        vec1 = Vector(*(p1-p0).data())
        vec2 = Vector(*(p2-p0).data())
        rot_axis = cross(vec1, vec2)
        ang = angle(vec1, vec2)/2 - 180
        rotation = Xform.rotation(rot_axis, ang)
        vec = rotation.apply(vec1)
        return vec

    return Vector(0, 0, 0)
    # Add more cases...?

    # else:
    #     # Returns this null Vector to ensure that the sinus perfomed below gets
    #     # equal to 0.0. This means that there is no preferential bonding
    #     # direction for this ligand
    #     return None

class LigandHelper(object):
    """
    Class that defines useful information about the metal's surrounding
    ligands.

    ligand : chimera.Atom
        Ligand being analyzed
    ligands : list of chimera.Atom
        All ligands in the surroundings, including `ligand`, 
        needed to test if `ligand` is bidentate or not.
    metal_center : chimera.Atom
        The metal ion subject to geometry analysis.
    """
    def __init__(self, ligand, ligands, metal_center):
        self.ligand = ligand
        self.metal_center = metal_center
        self.coord = self.bond_dir_origin = ligand.xformCoord()
        self.vector = self.coord - metal_center.xformCoord()
        self.neighbors = ligand.neighbors
        self.terminal = len(self.neighbors) == 1
        self.bidentate = False
        self.bidentate_mate = None

        for potential_mate in ligands:
            conditions = (potential_mate is not ligand, 
                          len(ligand.neighbors) == 1,
                          len(potential_mate.neighbors) == 1, 
                          ligand.neighbors[0] is potential_mate.neighbors[0])

            if all(conditions):
                self.bidentate = True
                self.bond_dir_origin = self.neighbors[0].xformCoord()
                self.bidentate_mate = self.potential_mate
                break