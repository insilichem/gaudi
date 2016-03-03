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
This module allows to explore molecular folding through normal modes analysis.

It works by calculating normal modes for the input molecule and moving along
a comgination of normal modes.

It needs at least a :class:`gaudi.genes.rotamers.molecule.Molecule`.

"""

# Python
import random
import logging
import numpy
# Chimera
import chimera
# Prody
import prody
# GAUDI
from gaudi.genes import GeneProvider
from gaudi import parse


ZERO = chimera.Point(0.0, 0.0, 0.0)
IDENTITY = ((1.0, 0.0, 0.0, 0.0),
            (0.0, 1.0, 0.0, 0.0),
            (0.0, 0.0, 1.0, 0.0))
logger = logging.getLogger(__name__)


def enable(**kwargs):
    # kwargs = NormalModes.validate(kwargs)
    return NormalModes(**kwargs)


class NormalModes(GeneProvider):

    """
    Parameters
    ----------
    molecule : chimera.molecule

    algorithm : callable, optional, default=None
        coarseGrain(prm) wich make mol.select().setBetas(i) where i
        is the index Coarse Grain group
        Where prm is prody AtomGroup

    alg_options : dictionary, optional
        Expected
        residues_number : int, optional, default=7
            number of residues per group
        mass_division : int, optional, default=100
            number of groups

    samples : int, optional, default=10000
        number of conformations to generate

    rmsd : float, optional, default=1.0
        average RMSD that the conformations will have with respect to the initial conformation


    Attributes
    ----------
    NORMAL_MODES : prody.modes
        normal modes calculated for the molecle
        stored in a prody modes class (ANM or RTB)

    NORMAL_MODE_SAMPLES : prody.ensemble
        configurations applying modes to molecule

    _original_coords : numpy.array
        Parent coordinates

    Notes
    -----

    """

    # validate = parse.Schema({
    #     'target': parse.Molecule()
    # }, extra=parse.ALLOW_EXTRA)

    NORMAL_MODES = None
    NORMAL_MODE_SAMPLES = None

    def __init__(self, molecule=None, algorithm=None, alg_options=None,
                 samples=10000, rmsd=1.0,
                 **kwargs):
        self.molecule = molecule
        self._original_coords = coords2numpy(molecule)
        self.algorithm = algorithm
        self.samples = samples
        self.rmsd = rmsd
        if alg_options:
            self._alg_options = alg_options
        else:
            self._alg_options = {}
        self._chimera_prody_dictionary = {}
        GeneProvider.__init__(self, **kwargs)

    def __ready__(self):
        if not self.NORMAL_MODES:
            self.prody_normal_modes()
        # self.allele = self.mutate(1.0)
        self.allele = random.choice(self.NORMAL_MODE_SAMPLES)

    def express(self):
        """
        aplicar cambio de coord
        """
        e = self._chimera_prody_dictionary
        for atom in self.molecule.atoms:
            new_coords = self.allele.getCoords()[e[atom.coordIndex]]
            atom.setCoord(chimera.Point(new_coords[0], new_coords[1], new_coords[2]))

    def unexpress(self):
        """
        revertir cambio coord
        """
        for i, atom in enumerate(self.molecule.atoms):
            original_coords = self._original_coords[i]
            atom.setCoord(
                chimera.Point(original_coords[0], original_coords[1], original_coords[2]))

    def mate(self, mate):
        """
        Combine coords between two samples in NORMAL_MODES_SAMPLES?
                            Or two samples between diferent NORMAL_MODES_SAMPLES?
        Or combine samples between two NORMAL_MODES_SAMPLES?

        De moment : pass
        """
        pass

    def mutate(self, indpb):
        if random.random() < self.indpb:
            return random.choice(self.NORMAL_MODE_SAMPLES)

    #####
    def prody_normal_modes(self):
        self.NORMAL_MODES, self._chimera_prody_dictionary, self._moldy = calc_normal_modes(
            self.molecule, self.algorithm, **self._alg_options)
        self.NORMAL_MODE_SAMPLES = prody.sampleModes(modes=self.NORMAL_MODES,
                                                     atoms=self._moldy,
                                                     n_confs=self.samples,
                                                     rmsd=self.rmsd)


####
def calc_normal_modes(mol, algorithm=None, **options):
    """
    Parameters
    ----------
    mol : chimera.Molecule or prody.AtomGroup
    algorithm : callable, optional, default=None
        coarseGrain(prm) wich make mol.select().setBetas(i) where i
        is the index Coarse Grain group
        Where prm is prody AtomGroup
    options : dict, optional
        Parameters for algorithm callable

    Returns
    -------
    modes : ProDy modes ANM or RTB
    e : dict
        dictionary: e[chimera_atom.coordIndex] = i-thm element prody getCoords() array
    moldy : prody.AtomGroup()
        ProDy molecule. Pass it to prody.sampleModes()
    """
    e = None
    if isinstance(mol, chimera.Molecule):
        moldy, e = chimera2prody(mol)
    elif isinstance(mol, prody.AtomGroup):
        moldy = mol

    modes = None
    if algorithm is not None:
        moldy = algorithm(moldy, **options)
        modes = prody.RTB('normal modes for {} using algorithm {}'.format(moldy.getTitle(), 1))
        modes.buildHessian(moldy.getCoords(), moldy.getBetas())
        modes.calcModes()
    else:
        modes = prody.ANM('normal modes for {}'.format(moldy.getTitle()))
        modes.buildHessian(moldy)
        modes.calcModes()
    return modes, e, moldy


def chimera2prody(mol):
    """
    Function that transforms a chimera molecule into a prody atom group

    Parameters
    ----------
    mol : chimera.Molecule

    Returns
    -------
    moldy : prody.AtomGroup()
    e : dict
        dictionary: e[chimera_atom.coordIndex] = i-thm element prody getCoords() array
    """
    moldy = prody.AtomGroup()
    try:
        coords, elements, names, resnums, chids, betas, masses = [], [], [], [], [], [], []
        d, e = {}, {}

        for i, atm in enumerate(mol.atoms):
            d[i] = atm.coordIndex
            e[atm.coordIndex] = i

            coords.append(tuple(atm.coord()))  # array documentation to improve
            elements.append(atm.element.name)
            names.append(atm.name)
            resnums.append(atm.residue.id.position)
            chids.append(atm.residue.id.chainId)
            masses.append(atm.element.mass)
            betas.append(atm.bfactor)

        moldy.setCoords(coords)
        moldy.setElements(elements)
        moldy.setNames(names)
        moldy.setResnums(resnums)
        moldy.setChids(chids)
        moldy.setBetas(betas)
        moldy.setMasses(masses)
        moldy.setTitle(str(mol.name))

        moldy.setBonds([[e[bond.atoms[0].coordIndex], e[bond.atoms[1].coordIndex]]
                        for bond in mol.bonds])

        moldy.setTitle(mol.name)

    except AttributeError:
        raise AttributeError('mol must be a chimera.Molecule')
    return moldy, e


def alg1(moldy, residues_number=7, **kwargs):
    """
    Coarse Grain Algorithm 1: groups per residues

    Parameters
    ----------
    moldy : prody.AtomGroup
    residues_number : int, optional, default=7
        number of residues per group

    Returns
    ------
    moldy : prody.AtomGroup
        New betas added
    """
    n = residues_number
    group = 1
    for chain in moldy.iterChains():
        num_residues = sorted(list(set(chain.getResnums())))
        chain_name = chain.getChid()
        for a, b in chunker(len(num_residues), n):
            try:
                start, end = num_residues[a-1], num_residues[b-1]
                selector = 'chain {} and resnum {} to {}'.format(chain_name, start, end)
                selection = moldy.select(selector)
                selection.setBetas(group)
                group += 1
            except AttributeError as e:
                print('Warning: {}'.format(e))
                pass
    return moldy


def alg2(moldy, mass_division=100, **kwargs):
    """
    Coarse Grain Algorithm 2: groups per mass percentage

    Parameters
    ----------
    moldy : prody.AtomGroup
    mass_division : int, optional, default=100
        number of groups

    Returns
    -------
    moldy: prody.AtomGroup
        New Betas added
    """
    group = 1

    M = sum(moldy.getMasses())
    m = M/mass_division
    mass = None

    for chain in moldy.iterChains():
        selection = moldy.select('chain {}'.format(chain.getChid()))
        mass = 0.

        for atom in iter(selection):
            atom.setBeta(group)
            mass += atom.getMass()
            if mass > m:
                mass = 0.
                group += 1
        group += 1
    return moldy


def alg3(moldy, max_bonds=3, **kwargs):
    group = 1

    for chain in moldy.iterChains():
        selection = moldy.select('chain {}'.format(chain.getChid()))
        for atom in iter(selection):
            atom.setBeta(group)
            if atom.numBonds() >= max_bonds:
                group += 1
        group += 1
    return moldy


def chunker(end, n):
    for i in range(0, end-n+1, n):
        yield i+1, i+n
    if end % n:
        yield end-end % n+1, end


def coords2numpy(molecule):
    """
    Parameters
    ----------
    molecule : chimera.molecule

    Returns
    -------
    numpy.array with molecule.atoms coordinates
    """
    return numpy.array([tuple(atom.coord()) for atom in molecule.atoms], dtype=float)
