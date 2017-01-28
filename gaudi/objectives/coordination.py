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
import numpy as np
# Chimera
import chimera
# GAUDI
from gaudi.objectives import ObjectiveProvider
from gaudi.exceptions import AtomsNotFound, ResiduesNotFound, MoleculesNotFound
from gaudi import parse
from gaudi._cpdrift import coherent_point_drift


logger = logging.getLogger(__name__)


GEOMETRIES = {
    'cube': np.array([(1.0, 1.0, 1.0), (1.0, 1.0, -1.0), (1.0, -1.0, 1.0), (1.0, -1.0, -1.0), (-1.0, 1.0, 1.0), (-1.0, 1.0, -1.0), (-1.0, -1.0, 1.0), (-1.0, -1.0, -1.0), (0.0, 0.0, 0.0)]),
    'hexagonal bipyramid': np.array([(0.0, 1.0, 0.0), (0.0, -1.0, 0.0), (0.0, 0.0, 1.0), (0.0, 0.0, -1.0), (0.8660, 0.0, 0.5), (0.8660, 0.0, -0.5), (-0.8660, 0.0, 0.5), (-0.8660, 0.0, -0.5), (0.0, 0.0, 0.0)]),
    'linear': np.array([(0.0, 1.0, 0.0), (0.0, -1.0, 0.0), (0.0, 0.0, 0.0)]),
    'octahedron': np.array([(1.0, 0.0, 0.0), (-1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, -1.0, 0.0), (0.0, 0.0, 1.0), (0.0, 0.0, -1.0), (0.0, 0.0, 0.0)]),
    'pentagonal bipyramid': np.array([(0.0, 1.0, 0.0), (0.0, -1.0, 0.0), (0.0, 0.0, 1.0), (0.9511, 0.0, 0.30901), (0.5878, 0.0, -0.8090), (-0.5878, 0.0, -0.8090), (-0.9511, 0.0, 0.3090), (0.0, 0.0, 0.0)]),
    'square planar': np.array([(1.0, 0.0, 0.0), (-1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, -1.0, 0.0), (0.0, 0.0, 0.0)]),
    'square pyramid': np.array([(0.0, 1.0, 0.0), (-0.6123, -0.5, 0.6123), (0.6123, -0.5, 0.6123), (-0.6123, -0.5, -0.6123), (0.6123, -0.5, -0.6123), (0.0, 0.0, 0.0)]),
    'tetrahedral': np.array([(1.0, 1.0, 1.0), (-1.0, -1.0, 1.0), (-1.0, 1.0, -1.0), (1.0, -1.0, -1.0), (0.0, 0.0, 0.0)]),
    'trigonal planar': np.array([(1.0, 0.0, 0.0), (-0.5, 0.0, 0.8660), (-0.5, 0.0,-0.8660), (0.0, 0.0, 0.0)]),
    'trigonal bipyramid': np.array([(0.0, 1.0, 0.0), (0.0, -1.0, 0.0), (1.0, 0.0, 0.0), (-0.5, 0.0, 0.8660), (-0.5, 0.0, -0.8660), (0.0, 0.0, 0.0)]),
    'trigonal prism': np.array([(-0.6547, 0.6547, 0.3779), (0.6547, 0.6547, 0.3779), (0.0, 0.6547, -0.7559), (-0.6547, -0.6547, 0.3779), (0.6547, -0.6547, 0.3779), (0.0, -0.6547, -0.7559), (0.0, 0.0, 0.0)])}

def enable(**kwargs):
    kwargs = Coordination.validate(kwargs)
    return Coordination(**kwargs)


class Coordination(ObjectiveProvider):

    """
    Coordination class

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
    _validate = {
        parse.Required('probe'): parse.Named_spec("molecule", "atom"),
        'radius': parse.Coerce(float),
        'atom_types': [basestring],
        'atom_names': [basestring],
        'atom_elements': [basestring],
        'residues': [parse.Named_spec("molecule", "residue")],
        'distance': parse.All(parse.Coerce(float), parse.Range(min=0)),
        'angle': parse.Coerce(float),
        'min_atoms': parse.All(parse.Coerce(int), parse.Range(min=2)),
        'geometry': parse.Any(parse.In(GEOMETRIES.keys()), [parse.Coordinates]),
        'enforce_all_residues': parse.Coerce(bool),
        'only_one_ligand_per_residue': parse.Coerce(bool),
        'prevent_intruders': parse.Coerce(bool),
        'method': parse.In(['simple', 'metalgeom', 'metalgeom_directional', 'cpd'])
        }
    
    def __init__(self, method='metalgeom_directional', probe=None, radius=None, atom_types=(),
                 atom_elements=(), atom_names=(), residues=(), geometry='tetrahedral',
                 distance=0, angle=None, dihedral=None, min_atoms=1, prevent_intruders=True,
                 enforce_all_residues=False, only_one_ligand_per_residue=False,
                 *args, **kwargs):
        ObjectiveProvider.__init__(self, **kwargs)
        self.method = method
        self._probe = probe
        self._residues = residues
        self.radius = radius
        self.atom_types = atom_types
        self.atom_names = atom_names
        self.atom_elements = atom_elements
        self.distance = distance
        self.angle = angle
        self.dihedral = dihedral
        self.min_atoms = min_atoms
        self.only_one_ligand_per_residue = only_one_ligand_per_residue
        self.enforce_all_residues = enforce_all_residues
        self.prevent_intruders = prevent_intruders
        if isinstance(geometry, basestring):
            self.geometry = np.copy(GEOMETRIES[geometry])
        else:
            self.geometry = np.array(geometry)
        self.n_vertices = self.geometry.shape[0] - 1
        if self.n_vertices < self.min_atoms:
            self.min_ligands = self.n_vertices
            logger.warn('# Vertices in selected geometry < min_ligands! Overriding '
                        'min_ligands with {}'.format(self.n_vertices))

    def molecules(self, ind):
        return [m.compound.mol for m in ind._molecules.values()]

    def probe(self, ind):
        mol, serial = self._probe
        return ind.find_molecule(mol).find_atom(serial)

    def residues(self, ind):
        for mol, pos in self._residues:
            for residue in ind.find_molecule(mol).find_residues(pos):
                yield residue

    def evaluate(self, ind):
        """
        1. Get requested atoms sorted by distance
        2. If they meet the minimum quantity, return the rmsd for
           that geometry
        3. If that's not possible of they are not enough, return penalty
        """
        try:
            test_atoms = [a for d, a in self.coordination_sphere(ind)]
        except ResiduesNotFound:
            logger.warning("Not enough atoms or some residues missing")
            return -1000 * self.weight
       
        missing_atoms = self.min_atoms - len(test_atoms)
        if missing_atoms > 0:
            logger.warning("Could not find enough ligand atoms in probe environment. "
                           "{} missing".format(missing_atoms))
            return -100 * missing_atoms * self.weight

        geometry = np.copy(self.geometry)
        ligands = test_atoms[:self.n_vertices]
        metal = self.probe(ind)
        atom_points = np.array([a.xformCoord() for a in ligands + [metal]])
        # rmsd
        try:
            _, _, rmsd = coherent_point_drift(atom_points, geometry, method='rigid',
                                              guess_steps=2, max_iterations=10)
        except Exception as e:
            logger.exception(e)  #
            logger.warning("Geometry not feasible in current conditions")
            return -1000 * self.weight
        # directionality
        directionality = sum(ideal_bond_deviation(metal, ligand, ligands) for ligand in ligands)
        
        return rmsd + directionality

    def coordination_sphere(self, ind):
        """
        1. Get atoms and residues found within `self.radius` angstroms from `self.probe`
        1.1. Found residues MUST include self.residues. Otherwise, apply penalty
        2. Sort atoms by absolute difference of `self.distance` and distance to `self.probe`.
           That way, nearest atoms are computed first.
        2.1. If found atoms do not include some of the requested types, apply penalty.
        """
        # Helpers
        def abs_distance(a):
            return abs(self.distance - metal.xformCoord().distance(a.xformCoord()))
        if self.atom_types:
            def atom_is_valid(a): return a.idatmType in self.atom_types
        elif self.atom_elements:
            def atom_is_valid(a): return a.element.name in self.atom_elements
        elif self.atom_names:
            def atom_is_valid(a): return a.name in self.atom_names
        else:
            def atom_is_valid(a): return True

        self._update_zone(ind)
        metal = self.probe(ind)
        residues = list(self.residues(ind))
        atoms = [a for a in self.zone.atoms() if a is not metal]
            
        atoms_by_distance = []
        found_residues = set()
        distance_and_atoms = sorted((abs_distance(a), a) for a in atoms)
        for d, a in distance_and_atoms:
            if atom_is_valid(a) and a.residue in residues and d > 1.0:
                atoms_by_distance.append((d, a))
                found_residues.add(a.residue)
            elif self.prevent_intruders:
                break

        if self.enforce_all_residues and found_residues != residues:
            logger.warning("Some atoms found, but some residues are missing")
            raise ResiduesNotFound

        return atoms_by_distance

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


def ideal_bond_deviation(metal, ligand, other_ligands=()):

    """
    Assess if the current bond vector is well oriented with
    respect to the ideal bond vector.

    Parameters
    ----------
    metal : chimera.Atom
        The ion `ligands` are coordinating to
    ligand : chimera.Atom
        Potential ligand atoms to `metal`

    Returns
    -------
    float
        Absolute sine of the angle between the ideal vector and
        the ligand-metal one.
    """
    ligand_idatm = chimera.idatm.typeInfo.get(ligand.idatmType)
    ligand_geometry = ligand_idatm.geometry if ligand_idatm else 3
    ideal_positions = ideal_bonded_positions(ligand, metal.element, geometry=ligand_geometry)
    if not ideal_positions:
        logger.warning("Ligand %s reports no available bonding positions. Check your atom_types!", ligand)
        return 1.0
    # else we can go on
    # Conditions and booleans
    rotates = ligand_geometry == 4 and len(ligand.neighbors) == 1
    bidentate = False
    for bidentate_mate in other_ligands:
        if bidentate_mate is not ligand:
            shared_neighbors = set(ligand.neighbors) & set(bidentate_mate.neighbors)
            if shared_neighbors:
                bidentate = True
                break

    # Coordinates and positions
    ligand_coord = ligand.xformCoord()
    metal_coord = metal.xformCoord()
    if not ligand.neighbors:
        # If ligand has no neighbors, it's an isolated atom. This means it will always 
        # be well oriented towards the metal. Then we can just return 0.0
        return 0.0
    # else we can go on
    neighbor = ligand.neighbors[0]
    neighbor_coord = neighbor.xformCoord()
    try:
        n_neighbor = next(a for a in neighbor.neighbors if a is not ligand)
        n_neighbor_coord = n_neighbor.xformCoord()
    except StopIteration:
        logger.warning('Ligand %s has no 2nd level neighbors. Dihedrals cannot be evaluated', ligand)
        rotates = True # This will force to ignore dihedral calculation ;)
    
    ideal_pos = ideal_positions[0]
    actual_angle = chimera.angle(neighbor_coord, ligand_coord, metal_coord)
    ideal_angles = [chimera.angle(neighbor_coord, ligand_coord, ideal_pos)]
    ideal_dihedrals = []
    if bidentate and not rotates:
        a = ligand_coord
        b = bidentate_mate.xformCoord() 
        c = shared_neighbors.pop().xformCoord()
        bidentate_ideal_pos = c + 1.5 * ((a-c) + (b-c))
        bidentate_ideal_angle = chimera.angle(neighbor_coord, ligand_coord, bidentate_ideal_pos)
        bidentate_ideal_dihedral = chimera.dihedral(n_neighbor_coord, neighbor_coord, 
                                                    ligand_coord, bidentate_ideal_pos)
        ideal_angles.append(bidentate_ideal_angle)
        ideal_dihedrals.append(bidentate_ideal_dihedral)
        
    angle_diff = min(delta - actual_angle for delta in ideal_angles)
    abs_sin_angle = abs(math.sin(math.radians(angle_diff)))
    if rotates:  # we don't care about the dihedral in this case
        return abs_sin_angle
    # else:
    actual_dihedral = chimera.dihedral(n_neighbor_coord, neighbor_coord, ligand_coord, metal_coord)
    ideal_dihedral = chimera.dihedral(n_neighbor_coord, neighbor_coord, ligand_coord, ideal_pos)
    ideal_dihedrals.append(ideal_dihedral)
    dihedral_diff = min(delta - actual_dihedral for delta in ideal_dihedrals)
    abs_sin_dihedral = abs(math.sin(math.radians(dihedral_diff)))
    return abs_sin_angle + abs_sin_dihedral


def ideal_bonded_positions(atom, element, geometry=None):
    if geometry is None:
        try:
            geometry = chimera.idatm.typeInfo[atom.idatmType].geometry
        except KeyError:
            geometry = 3

    bond_length = chimera.Element.bondLength(atom.element, element)
    neighbors_crd = [a.xformCoord() for a in atom.neighbors]
    return chimera.bondGeom.bondPositions(atom.xformCoord(), geometry, bond_length, neighbors_crd)
