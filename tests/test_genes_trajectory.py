#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
import mdtraj
import numpy as np
from conftest import datapath, expressed
from gaudi.genes.molecule import Molecule
from gaudi.genes.trajectory import Trajectory


def trajectory_setup(individual, molecule, traj, frame, max_frame):
    individual.genes['Molecule'] = Molecule(parent=individual, path=datapath(molecule))
    gene = Trajectory(parent=individual, path=datapath(traj), max_frame=max_frame,
                      target='Molecule')
    individual.genes['Trajectory'] = gene
    individual.__ready__()
    individual.__expression_hooks__()
    gene.allele = frame
    return gene


@pytest.mark.parametrize("mol_path, trajectory_path, frame, max_frame", [
    ('1amb.pdb', '1amb.dcd', 0, 1),
])
def test_trajectory(individual, mol_path, trajectory_path, frame, max_frame):
    gene = trajectory_setup(individual, mol_path, trajectory_path, frame, max_frame)
    traj = mdtraj.load_frame(datapath(trajectory_path), frame, top=datapath(mol_path))
    m = gene.molecule
    with expressed(individual):
        assert np.abs((m.compound.mol.activeCoordSet.xyzArray() - m.xyz()) < 0.001).all()
        assert np.abs((m.compound.mol.atomCoordinatesArray() - m.xyz()) < 0.001).all()
        assert np.abs((m.xyz() - (traj.xyz[0] * 10)) < 0.001).all()
        assert gene.topology.n_atoms == traj.top.n_atoms
        assert gene.topology.n_bonds == traj.top.n_bonds
        assert gene.topology.n_residues == traj.top.n_residues


@pytest.mark.parametrize("mol_path, trajectory_path, frame, max_frame", [
    ('1amb.pdb', '1amb.dcd', 0, 1),
])
def test_benchmark_rotamers(benchmark, individual, mol_path, trajectory_path,
                            frame, max_frame):
    @benchmark
    def run():
        trajectory_setup(individual, mol_path, trajectory_path, frame, max_frame)
