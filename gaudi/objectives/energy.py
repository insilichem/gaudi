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
This gene is a wrapper around OpenMM, providing a
GPU-accelerated energy calculation of the system
with a tiny molecular dynamics essay.
"""

# Python
try:
    from cStringIO import StringIO
except ImportError:
    from StringIO import StringIO
# 3rd party
import chimera
import simtk.openmm.app as openmm_app
from simtk import unit, openmm
# GAUDI
from gaudi.objectives import ObjectiveProvider


def enable(**kwargs):
    return Energy(**kwargs)


class Energy(ObjectiveProvider):

    """
    Calculate the energy of a system

    Parameters
    ----------
    forcefields : list of str
        Which forcefields to use

    Returns
    -------
    energy : float
        The estimated potential energy, in kJ/mol

    """

    def __init__(self, forcefields, *args, **kwargs):
        ObjectiveProvider.__init__(self, **kwargs)
        self.forcefields = forcefields
        pass

    def evaluate(self, individual):
        molecules = self.molecules(individual)
        topology, coordinates = self._chimera_molecule_to_openmm(molecules)
        energy = self._openmm_calculate_energy(topology, coordinates)

        return energy

    def molecules(self, individual):
        return tuple(m.compound.mol for m in individual.genes.values()
                     if m.__class__.__name__ == "Molecule")

    def _openmm_calculate_energy(self, topology, coordinates):
        """
        Set up an OpenMM simulation with default parameters
        and return the potential energy of the initial state

        Parameters
        ----------
        topology : simtk.openmm.app.topology.Topology
            Topology of the system, in an OpenMM compatible format
        coordinates : simtk.unit.Quantity
            Positions of the atoms in the system

        Returns
        -------
        potential_energy : float
            Potential energy of the system, in kJ/mol
        """
        # Prepare system
        forcefield = openmm_app.ForceField(*self.forcefields)
        system = forcefield.createSystem(topology, nonbondedMethod=openmm_app.PME,
                                         nonbondedCutoff=1.0*unit.nanometers,
                                         constraints=openmm_app.HBonds, rigidWater=True,
                                         ewaldErrorTolerance=0.0005)
        integrator = openmm.LangevinIntegrator(300*unit.kelvin, 1.0/unit.picoseconds,
                                               2.0*unit.femtoseconds)
        integrator.setConstraintTolerance(0.00001)
        # Prepare simulation
        simulation = openmm_app.Simulation(topology, system, integrator)
        simulation.context.setPositions(coordinates)
        # Retrieve initial energy
        state = simulation.context.getState(getEnergy=True)
        potential_energy = state.getPotentialEnergy()

        return potential_energy

    @staticmethod
    def _chimera_molecule_to_openmm(molecules):
        """
        Convert a Chimera Molecule object to OpenMM structure,
        providing topology and coordinates.

        Parameters
        ----------
        molecule : chimera.Molecule

        Returns
        -------
        topology : simtk.openmm.app.topology.Topology
        coordinates : simtk.unit.Quantity

        """
        memfile = StringIO()
        chimera.pdbWrite(molecules, molecules[0].openState.xform, memfile)
        pdb = openmm_app.PDBFile(memfile)
        memfile.close()

        return pdb.getTopology(), pdb.getPositions()
