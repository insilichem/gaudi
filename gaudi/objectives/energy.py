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
This objective is a wrapper around OpenMM, providing a GPU-accelerated energy 
calculation of the system with a simple forcefield evaluation.
"""

# Python
try:
    from cStringIO import StringIO
except ImportError:
    from StringIO import StringIO
import os
# 3rd party
import chimera
import simtk.openmm.app as openmm_app
from simtk import unit, openmm
from openmoltools.amber import run_antechamber
from openmoltools.utils import create_ffxml_file
# GAUDI
from gaudi import parse
from gaudi.objectives import ObjectiveProvider

_openmm_builtin_forcefields = os.listdir(os.path.join(openmm_app.__path__[0], 'data'))

def enable(**kwargs):
    kwargs = Energy.validate(kwargs)
    return Energy(**kwargs)


class Energy(ObjectiveProvider):

    """
    Calculate the energy of a system

    Parameters
    ----------
    forcefields : list of str
        Which forcefields to use
    auto_parametrize: list of str
        List of Molecule instances GAUDI should try to auto parametrize with antechamber.

    Returns
    -------
    energy : float
        The estimated potential energy, in kJ/mol

    """

    validate = parse.Schema({
        'forcefields': parse.Any(parse.ExpandUserPathExists, parse.In(_openmm_builtin_forcefields)),
        'auto_parametrize': [parse.Molecule_name]
        }, extra=parse.ALLOW_EXTRA)

    def __init__(self, forcefields=('amber99sbildn.xml',), auto_parametrize=None, *args, **kwargs):
        ObjectiveProvider.__init__(self, **kwargs)
        self.auto_parametrize = auto_parametrize
        self.topology = None
        self._simulation = None

        if auto_parametrize:
            filenames = [g.path for m in auto_parametrize
                         for g in self.environment.cfg.genes
                         if g.name == m]
            additional_ffxml = self._gaff2xml(*filenames)
            forcefields.append(additional_ffxml)

        self.forcefields = forcefields
        self.forcefield = openmm_app.ForceField(*self.forcefields)

    def evaluate(self, individual):
        """
        Calculates the energy of current individual
        """
        molecules = self.molecules(individual)
        topology, coordinates = self.chimera_molecule_to_openmm(*molecules)
        # Topology changed -> rebuild universe
        # Instead of checking if the topologies are equivalent,
        # check if GAUDI input can produce more than one topology
        if not self._gaudi_is_static(individual):
            # This forces a Simulation rebuild
            self._simulation = None

        self.topology = topology
        energy = self.calculate_energy(coordinates)

        return energy

    def molecules(self, individual):
        return tuple(m.compound.mol for m in individual.genes.values()
                     if m.__class__.__name__ == "Molecule")

    @property
    def simulation(self):
        """
        Build a new OpenMM simulation if not yet defined and return it

        Notes
        -----
        self.topology must be defined previously! Use self.chimera_molecule_to_openmm()
        to get topology and coordinates, and set it manually.

        """
        if self._simulation is None:
            system = self.forcefield.createSystem(self.topology,
                                                  nonbondedMethod=openmm_app.CutoffNonPeriodic,
                                                  nonbondedCutoff=1.0*unit.nanometers,
                                                  rigidWater=True, constraints=None)
            integrator = openmm.VerletIntegrator(0.001)
            self._simulation = openmm_app.Simulation(self.topology, system, integrator)
        return self._simulation

    def calculate_energy(self, coordinates):
        """
        Set up an OpenMM simulation with default parameters
        and return the potential energy of the initial state

        Parameters
        ----------
        coordinates : simtk.unit.Quantity
            Positions of the atoms in the system

        Returns
        -------
        potential_energy : float
            Potential energy of the system, in kJ/mol
        """
        self.simulation.context.setPositions(coordinates)
        # Retrieve initial energy
        state = self.simulation.context.getState(getEnergy=True)
        potential_energy = state.getPotentialEnergy()

        return potential_energy._value

    @staticmethod
    def chimera_molecule_to_openmm(*molecules):
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
        # Write PDB to StringIO (memory file)
        pdbfile = StringIO()
        for molecule in molecules:
            chimera.pdbWrite([molecule], molecule.openState.xform, pdbfile)
        pdbfile.seek(0)

        molecule = openmm_app.PDBFile(pdbfile)
        return molecule.topology, molecule.positions

    @staticmethod
    def _gaff2xml(*filenames):
        """
        Use OpenMolTools wrapper to run antechamber programatically
        and auto parametrize requested molecules.

        Parameters
        ----------
        filenames: list of str
            List of the filenames of the molecules to parametrize

        Returns
        -------
        ffxmls : StringIO
            Compiled ffxml file produced by antechamber and openmoltools converter
        """
        frcmods, gaffmol2s = [], []
        for filename in filenames:
            name = '.'.join(filename.split('.')[:-1])
            gaffmol2, frcmod = run_antechamber(name, filename)
            frcmods.append(frcmod)
            gaffmol2s.append(gaffmol2)
        return create_ffxml_file(gaffmol2s, frcmods)


    def _gaudi_is_static(self, individual):
        """
        Check if this essay is performing topology changes.

        Genes that can change topologies:
            - gaudi.genes.rotamers with mutations ON
            - gaudi.genes.molecule with block building enabled

        Parameters
        ----------
        individual : gaudi.base.Individual
            The individual to be analyzed for dynamic behaviour

        Returns
        -------
        bool
        """
        for gene in individual.genes.values():
            if gene.__class__.__name__ == 'Rotamers':
                if gene.mutations:
                    return False
            if gene.__class__.__name__ == 'Molecule':
                if len(gene.catalog) > 1:
                    return False

        return True


def calculate_energy(filename, forcefield=None):
    """
    Calculate energy from PDB file with desired forcefield. If not specified,
    amber99sbildn will be used. Returns potential energy in kJ/mol.
    """
    if forcefield is None:
        forcefield = ['amber99sbildn.xml']

    mol = openmm_app.PDBFile(filename)
    ff = openmm_app.ForceField(*forcefield)
    system = ff.createSystem(mol.topology,
                             nonbondedMethod=openmm_app.CutoffNonPeriodic,
                             nonbondedCutoff=1.0*unit.nanometers,
                             rigidWater=True, constraints=None)
    integrator = openmm.VerletIntegrator(0.001)
    sim = openmm_app.Simulation(mol.topology, system, integrator)
    sim.context.setPositions(mol.positions)
    state = sim.context.getState(getEnergy=True)
    potential_energy = state.getPotentialEnergy()

    return potential_energy._value


def _test_topology_equality(t1, t2):
    """
    Test if two topologies are identical

    Parameters
    ----------
    t1, t2 : simtk.openmm.app.topology.Topology

    Returns
    -------
    bool

    Notes
    -----
    Ripped straight from ``mdtraj.core.topology.__eq__``. It might be slow, though.
    """
    if not isinstance(t1, openmm_app.topology.Topology):
        return False
    if not isinstance(t2, openmm_app.topology.Topology):
        return False

    if t1 is t2:
        return True

    if len(t1._chains) != len(t2._chains):
        return False

    for c1, c2 in zip(t1.chains(), t2.chains()):
        if c1.index != c2.index:
            return False
        if len(c1._residues) != len(c2._residues):
            return False

        for r1, r2 in zip(c1.residues(), c2.residues()):
            if (r1.index != r1.index) or (r1.name != r2.name):  # or (r1.resSeq != r2.resSeq):
                return False
            if len(r1._atoms) != len(r2._atoms):
                return False

            for a1, a2 in zip(r1.atoms(), r2.atoms()):
                if (a1.index != a2.index) or (a1.name != a2.name):
                    return False
                if a1.element is not None and a2.element is not None:
                    if a1.element != a2.element:
                        return False
                    # for attr in ['atomic_number', 'name', 'symbol']:
                    #    if getattr(a1.element, attr) != getattr(a2.element, attr):
                    #        return False

    if len(t1._bonds) != len(t2._bonds):
        return False

    # the bond ordering is somewhat ambiguous, so try and fix it for comparison
    t1_sorted_bonds = sorted([(a1.index, b1.index) for (a1, b1) in t1.bonds()])
    t2_sorted_bonds = sorted([(a2.index, b2.index) for (a2, b2) in t2.bonds()])

    for i in range(len(t1._bonds)):
        (a1, b1) = t1_sorted_bonds[i]
        (a2, b2) = t2_sorted_bonds[i]
        if (a1 != a2) or (b1 != b2):
            return False

    return True

if __name__ == "__main__":
    import sys
    print(calculate_energy(sys.argv[1]))