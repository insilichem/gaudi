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
This objective calculates SASA for the given system (or region).
"""

# Python
import logging
from math import ceil
# 3rd party
import numpy as np
# Chimera
import chimera
from _surface import surface_area, enclosed_volume
from _multiscale import get_atom_coordinates, bounding_box
from _gaussian import sphere_surface_distance
from _contour import surface as contour_surface
from Matrix import transform_points
# GAUDI
from gaudi import parse
from gaudi.objectives import ObjectiveProvider

logger = logging.getLogger(__name__)


def enable(**kwargs):
    kwargs = Solvation.validate(kwargs)
    return Solvation(**kwargs)


class Solvation(ObjectiveProvider):

    """
    Solvation class

    Parameters
    ----------
    targets : [str]
        Names of the molecule genes being analyzed
    threshold : float
        Optimize the difference to this value
    radius : float
        Max distance to search for neighbor atoms from targets.
    """

    _validate = {
        parse.Required('targets'): [parse.Molecule_name],
        'threshold': parse.All(parse.Coerce(float), parse.Range(min=0)),
        'radius': parse.All(parse.Coerce(float), parse.Range(min=0)),
        'method': parse.In(['volume', 'area'])
        }

    def __init__(self, targets=None, threshold=0.0, radius=5.0, method='area',
                 *args, **kwargs):
        ObjectiveProvider.__init__(self, **kwargs)
        self._targets = targets
        self.threshold = threshold
        self.radius = radius
        self.method = method
        if method == 'area':
            self.evaluate = self.evaluate_area
        else:
            self.evaluate = self.evaluate_volume

    def targets(self, ind):
        return [ind.find_molecule(target).compound.mol for target in self._targets]

    def molecules(self, ind):
        return tuple(m.compound.mol for m in ind._molecules.values())

    def surface(self, ind):
        atoms = self.zone_atoms(self.targets(ind), self.molecules(ind))
        return grid_sas_surface(atoms)
        
    def evaluate_area(self, ind):
        return abs(surface_area(*self.surface(ind)) - self.threshold)

    def evaluate_volume(self, ind):
        return abs(enclosed_volume(*self.surface(ind))[0] - self.threshold)

    def zone_atoms(self, probes, molecules):
        self.zone.clear()
        self.zone.add([a for probe in probes for a in probe.atoms])
        if self.radius:
            self.zone.merge(chimera.selection.REPLACE,
                            chimera.specifier.zone(self.zone, 'atom', None, 
                                                   self.radius, molecules))
        return self.zone.atoms()


def grid_sas_surface(atoms, probe_radius=1.4, grid_spacing=0.5):
    """
    Stripped from Chimera's Surface.gridsurf
    """
    xyz = get_atom_coordinates(atoms, transformed = False)
    radii = np.array([a.radius for a in atoms], np.float32)

    # Compute bounding box for atoms
    xyz_min, xyz_max = bounding_box(xyz)
    pad = 2*probe_radius + radii.max()
    origin = [x-pad for x in xyz_min]

    # Create 3d grid for computing distance map
    s = grid_spacing
    shape = [int(ceil((xyz_max[a] - xyz_min[a] + 2*pad) / s))
             for a in (2,1,0)]
    matrix = np.empty(shape, np.float32)
    max_index_range = 2
    matrix[:,:,:] = max_index_range

    # Transform centers and radii to grid index coordinates
    xyz_to_ijk_tf = ((1.0/s, 0, 0, -origin[0]/s),
                     (0, 1.0/s, 0, -origin[1]/s),
                     (0, 0, 1.0/s, -origin[2]/s))
    ijk = xyz.copy()
    transform_points(ijk, xyz_to_ijk_tf)
    probed_radii = radii.copy()
    probed_radii += probe_radius
    probed_radii /= s

    # Compute distance map from surface of spheres, positive outside.
    sphere_surface_distance(ijk, probed_radii, max_index_range, matrix)
    # Get the SAS surface as a contour surface of the distance map
    return contour_surface(matrix, 0, cap_faces=False, calculate_normals=False)

