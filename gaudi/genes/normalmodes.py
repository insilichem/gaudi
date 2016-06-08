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
a combination of normal modes.

It needs at least a :class:`gaudi.genes.rotamers.molecule.Molecule`.

"""

# Python
from __future__ import print_function, division
import random
import logging
import numpy
import os
# Chimera
import chimera
# 3rd party
import prody
from repoze.lru import LRUCache
from cclib.parser import Gaussian
# GAUDI
from gaudi.genes import GeneProvider
from gaudi import parse


logger = logging.getLogger(__name__)


def enable(**kwargs):
    kwargs = NormalModes.validate(kwargs)
    return NormalModes(**kwargs)


class NormalModes(GeneProvider):

    """
    Parameters
    ----------
    method : str
        Expected
            prody : calculate normal modes using prody algorithms
            gaussian : read normal modes from a gaussian output file
    target : str
        Name of the Gene containing the actual molecule
    modes : list, optional, default=range(12)
        Modes to be used to move the molecule
    group_by : string, optional, default=None
        group_by_* algorithm name
        group_by_* : callable, optional, default=None
            coarseGrain(prm) wich make mol.select().setBetas(i) where i
            is the index Coarse Grain group
            Where prm is prody AtomGroup
    group_lambda : dictionary, optional
        Expected
        residues_number : int, optional, default=7
            number of residues per group
        mass_division : int, optional, default=100
            number of groups
    path : str
        Gaussian or prody modes output path
        Obligatory if method=gaussian
    write_modes: boolean, optional
        write a molecule_modes.nmd file with the pordy modes
    n_samples : int, optional, default=10000
        number of conformations to generate
    rmsd : float, optional, default=1.0
        average RMSD that the conformations will have with respect to the initial conformation

    Attributes
    ----------
    NORMAL_MODES : prody.modes
        normal modes calculated for the molecule or readed
        from the gaussian frequencies output file stored
        in a prody modes class (ANM or RTB)
    NORMAL_MODE_SAMPLES : prody.ensemble
        configurations applying modes to molecule
    _original_coords : numpy.array
        Parent coordinates
    _chimera2prody : dict
        _chimera2prody[chimera_index] = prody_index
    """

    validate = parse.Schema({
        parse.Required('method'): parse.In(['prody', 'gaussian']),
        'path': parse.RelPathToInputFile(),
        'write_modes': parse.Boolean,
        parse.Required('target'): parse.Molecule_name,
        'group_by': parse.In(['residues', 'mass', 'calpha', '']),
        'group_lambda': parse.All(parse.Coerce(int), parse.Range(min=1)),
        'modes': [parse.All(parse.Coerce(int), parse.Range(min=0))],
        'n_samples': parse.All(parse.Coerce(int), parse.Range(min=1)),
        'rmsd': parse.All(parse.Coerce(float), parse.Range(min=0))
    }, extra=parse.ALLOW_EXTRA)

    def __init__(self, method='prody', target=None, modes=None, n_samples=10000, rmsd=1.0,
                 group_by=None, group_lambda=None,
                 path=None, write_modes=False, **kwargs):
        # Fire up!
        GeneProvider.__init__(self, **kwargs)
        self.method = method
        self.target = target
        self.modes = modes if modes is not None else range(12)
        self.max_modes = max(self.modes) + 1
        self.n_samples = n_samples
        self.rmsd = rmsd
        self.group_by = None
        self.group_by_options = None
        self.path = None
        self.write_modes = write_modes
        if method == 'prody':
            if path is None:
                self.normal_modes_function = self.calculate_prody_normal_modes
                self.group_by = group_by
                self.group_by_options = {} if group_lambda is None else {'n': group_lambda}
            else:
                self.path = path
                self.normal_modes_function = self.read_prody_normal_modes
        else:  # gaussian
            self.normal_modes_function = self.read_gaussian_normal_modes
            if path is None:
                raise ValueError('Path is required if method == gaussian')
            self.path = path

        if self.name not in self._cache:
            self._cache[self.name] = LRUCache(300)

    def __ready__(self):
        """
        Second stage of initialization

        It saves the parent coordinates, calculates the normal modes and initializes the allele
        """
        cached = self._CACHE.get('normal_modes')
        if not cached:
            normal_modes, normal_modes_samples, chimera2prody, prody_molecule = self.normal_modes_function()
            self._CACHE.put('normal_modes', normal_modes)
            self._CACHE.put('normal_modes_samples', normal_modes_samples)
            self._CACHE.put('chimera2prody', chimera2prody)
            self._CACHE.put('original_coords', chimeracoords2numpy(self.molecule))
            if self.write_modes:
                title = os.path.join(self.parent.cfg.output.path,
                                     '{}_modes.nmd'.format(self.molecule.name))
                prody.writeNMD(title, normal_modes, prody_molecule)
        self.allele = random.choice(self.NORMAL_MODES_SAMPLES)

    def express(self):
        """
        Apply new coords as provided by current normal mode
        """
        c2p = self._chimera2prody
        for atom in self.molecule.atoms:
            index = c2p[atom.serialNumber]
            new_coords = self.allele[index]
            atom.setCoord(chimera.Point(*new_coords))

    def unexpress(self):
        """
        Undo coordinates change
        """
        for i, atom in enumerate(self.molecule.atoms):
            atom.setCoord(chimera.Point(*self._original_coords[i]))

    def mate(self, mate):
        """
        Combine coords between two samples in NORMAL_MODES_SAMPLES?
                            Or two samples between diferent NORMAL_MODES_SAMPLES?
        Or combine samples between two NORMAL_MODES_SAMPLES?

        For now : pass
        """
        pass

    def mutate(self, indpb):
        """
        (mutate to/get) another SAMPLE with probability = indpb
        """
        if random.random() < self.indpb:
            return random.choice(self.NORMAL_MODES_SAMPLES)

    #####
    @property
    def molecule(self):
        return self.parent.genes[self.target].compound.mol

    @property
    def _CACHE(self):
        return self._cache.get(self.name)

    @property
    def NORMAL_MODES(self):
        return self._CACHE.get('normal_modes')

    @property
    def NORMAL_MODES_SAMPLES(self):
        return self._CACHE.get('normal_modes_samples')

    @property
    def _chimera2prody(self):
        return self._CACHE.get('chimera2prody')

    @property
    def _original_coords(self):
        return self._CACHE.get('original_coords')

    def calculate_prody_normal_modes(self):
        """
        calculate normal modes, creates a diccionary between chimera and prody indices
        and calculate n_confs number of configurations using this modes
        """
        prody_molecule, chimera2prody = convert_chimera_molecule_to_prody(self.molecule)
        modes = prody_modes(prody_molecule, self.max_modes, GROUPERS[self.group_by],
                            **self.group_by_options)
        samples = prody.sampleModes(modes=modes[self.modes], atoms=prody_molecule,
                                    n_confs=self.n_samples, rmsd=self.rmsd)
        samples.addCoordset(prody_molecule)
        samples_coords = [sample.getCoords() for sample in samples]
        return modes, samples_coords, chimera2prody, prody_molecule

    def read_prody_normal_modes(self):
        prody_molecule, chimera2prody = convert_chimera_molecule_to_prody(self.molecule)
        modes = prody.parseNMD(self.path)[0]
        samples = prody.sampleModes(modes=modes[self.modes], atoms=prody_molecule,
                                    n_confs=self.n_samples, rmsd=self.rmsd)
        samples.addCoordset(prody_molecule)
        samples_coords = [sample.getCoords() for sample in samples]
        return modes, samples_coords, chimera2prody, prody_molecule

    def read_gaussian_normal_modes(self):
        """
        read normal modes, creates a diccionary between chimera and prody indices
        and calculate n_confs number of configurations using this modes
        """
        prody_molecule, chimera2prody = convert_chimera_molecule_to_prody(self.molecule)
        modes = gaussian_modes(self.path)

        samples = prody.sampleModes(modes=modes[self.modes], atoms=prody_molecule,
                                    n_confs=self.n_samples, rmsd=self.rmsd)
        samples.addCoordset(prody_molecule)
        samples_coords = [sample.getCoords() for sample in samples]
        return modes, samples_coords, chimera2prody, prody_molecule


####
def prody_modes(molecule, max_modes, algorithm=None, **options):
    """
    Parameters
    ----------
    molecule : prody.AtomGroup
    nax_modes : int
        number of modes to calculate
    algorithm : callable, optional, default=None
        coarseGrain(prm) wich make molecule.select().setBetas(i) where i
        is the index Coarse Grain group
        Where prm is prody AtomGroup
    options : dict, optional
        Parameters for algorithm callable

    Returns
    -------
    modes : ProDy modes ANM or RTB
    """
    modes = None
    if algorithm in ['residues', 'mass']:
        title = 'normal modes for {}'.format(molecule.getTitle())
        molecule = algorithm(molecule, **options)
        modes = prody.RTB(title)
        modes.buildHessian(molecule.getCoords(), molecule.getBetas())
        modes.calcModes(n_modes=max_modes)
    elif algorithm == 'calpha':
        calphas_modes = prody.ANM('normal modes for {}'.format(molecule.getTitle()))
        calphas = molecule = molecule.select(algorithm)
        calphas_modes.buildHessian(calphas)
        calphas_modes.calcModes(n_modes=max_modes)
        modes = prody.extendModel(calphas_modes, calphas, molecule, norm=True)[0]
    else:
        modes = prody.ANM('normal modes for {}'.format(molecule.getTitle()))
        modes.buildHessian(molecule)
        modes.calcModes(n_modes=max_modes)
    return modes


def gaussian_modes(path):
    """
    Read the modes
    Create a prody.modes instance

    Parameters
    ----------
    path : str
        gaussian frequencies output path

    Returns
    -------
    modes : ProDy modes ANM or RTB
    """
    gaussian_parser = Gaussian(path).parse()
    shape = gaussian_parser.vibdisps.shape
    modes_vectors = gaussian_parser.vibdisps.reshape(shape[0], shape[1]*shape[2]).T
    modes_frequencies = numpy.abs(gaussian_parser.vibfreqs)
    modes = prody.NMA()
    modes.setEigens(vectors=modes_vectors, values=modes_frequencies)
    return modes


def convert_chimera_molecule_to_prody(molecule):
    """
    Function that transforms a chimera molecule into a prody atom group

    Parameters
    ----------
    molecule : chimera.Molecule

    Returns
    -------
    prody_molecule : prody.AtomGroup()
    chimera2prody : dict
        dictionary: chimera2prody[chimera_atom.coordIndex] = i-thm element prody getCoords() array
    """
    prody_molecule = prody.AtomGroup()
    try:
        coords, elements, serials = [], [], []
        names, resnums, resnames = [], [], []
        chids, betas, masses = [], [], []
        chimera2prody = {}
        offset_chimera_residue = min(r.id.position for r in molecule.residues)

        for i, atm in enumerate(molecule.atoms):
            chimera2prody[atm.serialNumber] = i
            coords.append(tuple(atm.coord()))  # array documentation to improve
            elements.append(atm.element.name)
            serials.append(atm.serialNumber)
            names.append(atm.name)
            resnums.append(atm.residue.id.position - offset_chimera_residue)
            resnames.append(atm.residue.type)
            chids.append(atm.residue.id.chainId)
            masses.append(atm.element.mass)
            betas.append(atm.bfactor)

        prody_molecule.setCoords(coords)
        prody_molecule.setElements(elements)
        prody_molecule.setSerials(serials)
        prody_molecule.setNames(names)
        prody_molecule.setResnums(resnums)
        prody_molecule.setResnames(resnames)
        prody_molecule.setChids(chids)
        prody_molecule.setBetas(betas)
        prody_molecule.setMasses(masses)
        prody_molecule.setTitle(str(molecule.name))
        prody_molecule.setBonds([(chimera2prody[bond.atoms[0].serialNumber],
                                  chimera2prody[bond.atoms[1].serialNumber]) for bond in molecule.bonds])

    except AttributeError:
        raise TypeError('Attribute not found. Molecule must be a chimera.Molecule')

    return prody_molecule, chimera2prody


def group_by_residues(molecule, n=7):
    """
    Coarse Grain Algorithm 1: groups per residues

    Parameters
    ----------
    molecule : prody.AtomGroup
    n : int, optional, default=7
        number of residues per group

    Returns
    ------
    molecule : prody.AtomGroup
        New betas added
    """
    group = 1
    for chain in molecule.iterChains():
        residues_indices = sorted(list(set(chain.getResnums())))
        chain_name = chain.getChid()
        for a, b in chunker(len(residues_indices), n):
            try:
                start, end = residues_indices[a-1], residues_indices[b-1]
                selector = 'chain {} and resnum {} to {}'.format(chain_name, start, end)
                selection = molecule.select(selector)
                selection.setBetas(group)
                group += 1
            except AttributeError as e:
                logger.warning(str(e))
    return molecule


def group_by_mass(molecule, n=100):
    """
    Coarse Grain Algorithm 2: groups per mass percentage

    Parameters
    ----------
    molecule : prody.AtomGroup
    n: int, optional, default=100
        Intended number of groups. The mass of the system will be divided by this number,
        and each group will have the corresponding proportional mass. However, the final
        number of groups can be slightly different.

    Returns
    -------
    molecule: prody.AtomGroup
        New Betas added
    """
    group = 1

    total_mass = sum(molecule.getMasses())
    chunk_mass = total_mass/n

    for chain in molecule.iterChains():
        selection = molecule.select('chain {}'.format(chain.getChid()))
        mass_accumulator = 0.

        for atom in iter(selection):
            atom.setBeta(group)
            mass_accumulator += atom.getMass()
            if mass_accumulator > chunk_mass:
                mass_accumulator = 0.
                group += 1
        group += 1
    return molecule


def alg3(moldy, max_bonds=3, **kwargs):
    """
    TESTS PENDING!

    Coarse Grain Algorithm 3: Graph algorithm.
        New group when a vertice: have more than n,
                                  have 0 edges
                                  new chain

    Parameters
    ----------
    moldy : prody.AtomGroup
    n : int, optional, default=2
        maximum bonds number

    Returns
    -------
    moldy: prody.AtomGroup
        New Betas added
    """
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
    """
    divide end integers in closed groups of n
    """
    for i in range(0, end-n+1, n):
        yield i+1, i+n
    if end % n:
        yield end-end % n+1, end


def chimeracoords2numpy(molecule):
    """
    Parameters
    ----------
    molecule : chimera.molecule

    Returns
    -------
    numpy.array with molecule.atoms coordinates
    """
    return numpy.array([tuple(atom.coord()) for atom in molecule.atoms], dtype=float)


GROUPERS = {
    'residues': group_by_residues,
    'mass': group_by_mass,
    'calpha': 'calpha',
    '': None
}
