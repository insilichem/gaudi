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

It works by opening a gaussian file, converting vectors to prody modes and moving along
a combination of normal modes.

It needs at least a :class: `gaudi.genes.rotamers.molecule.Molecule`
                    :path: `path to gaussian file`

"""

# Python
from __future__ import print_function, division
import random
import logging
import numpy
# Chimera
import chimera
# 3rd party
import prody
from repoze.lru import LRUCache
# GAUDI
from gaudi.genes import GeneProvider
from gaudi import parse


logger = logging.getLogger(__name__)


def enable(**kwargs):
    kwargs = GaussianNormalModes.validate(kwargs)
    return GaussianNormalModes(**kwargs)


class GaussianNormalModes(GeneProvider):

    """
    Parameters
    ----------
    target : str
        Name of the Gene containing the actual molecule.

    n_modes : int, optional, default=5
        number of modes used to calculate conformations

    n_samples : int, optional, default=10000
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

    _chimera2prody : dict
        _chimera2prody[chimera_index] = prody_index
    
    Notes
    -----
    

    """

    validate = parse.Schema({
        'target': parse.Molecule_name,
        parse.Required('path'): parse.RelPathToInputFile(),
        # 'n_modes': parse.All(parse.Coerce(int), parse.Range(min=1)),
        'n_modes': [int],
        'n_samples': parse.All(parse.Coerce(int), parse.Range(min=1)),
        'rmsd': parse.All(parse.Coerce(float), parse.Range(min=0))
    }, extra=parse.ALLOW_EXTRA)

    def __init__(self, target=None, path=None, n_modes=range(6), n_samples=10000, rmsd=1.0,
                 **kwargs):
        # Fire up!
        GeneProvider.__init__(self, **kwargs)
        self.target = target
        self.path = path
        self.n_modes = n_modes
        self.n_samples = n_samples
        self.rmsd = rmsd
        if self.name not in self._cache:
            self._cache[self.name] = LRUCache(300)
        

    def __ready__(self):
        """
        Second stage of initialization

        It saves the parent coordinates, reads the file and initializes the allele
        """
        cached = self._CACHE.get('normal_modes')
        if not cached:
            normal_modes, normal_modes_samples, chimera2prody = self.calculate_normal_modes()
            self._CACHE.put('normal_modes', normal_modes)
            self._CACHE.put('normal_modes_samples', normal_modes_samples)
            self._CACHE.put('chimera2prody', chimera2prody)
            self._CACHE.put('original_coords', chimeracoords2numpy(self.molecule))
        # self.allele = self.mutate(1.0)
        self.allele = random.choice(self.NORMAL_MODES_SAMPLES)

    def express(self):
        """
        Apply new coords as provided by current normal mode
        """
        c2p = self._chimera2prody
        for atom in self.molecule.atoms:
            index = c2p[atom.coordIndex]
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
        # return self.target

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

    def calculate_normal_modes(self):
        """
        calculate normal modes, creates a diccionary between chimera and prody indices
        and calculate n_confs number of configurations using this modes
        """
        prody_molecule, chimera2prody = convert_chimera_molecule_to_prody(self.molecule)
        modes = gaussian_modes(self.path)
        samples = prody.sampleModes(modes=modes[self.n_modes], atoms=prody_molecule,
                                    n_confs=self.n_samples, rmsd=self.rmsd)
        samples_coords = [sample.getCoords() for sample in samples]
        return modes, samples_coords, chimera2prody


####
def gaussian_modes(path):
    """
    Read the modes and store as vectors
    Create a prody.modes instance with the vectors
    Check if the molecule in the modes is our molecule?
    """
    vector_modes = read_vectors(path)
    frequency_modes = read_frequencies(path)
    for i, freq in enumerate(frequency_modes):
        frequency_modes[i] = abs(freq)
    modes = create_prody_modes(frequency_modes, vector_modes)
    return modes


def read_vectors(path):
    """
    Read normal modes vectors from gaussian file

    Parameters
    ----------
    path : str

    Return
    ------
    modes : list

    Dubte: com sabem que en una linia hi trobarem tres trios de vectors?
    """
    input_file = open(path, "r")

    modes = []
    num_atoms = None
    # vfreq1 = []
    # vfreq2 = []
    # vfreq3 = []
    # vvoid = [float(0), float(0), float(0)]

    for line in input_file:

        if line.find('Input orientation') != -1 or line.find('Standard orientation') != -1 and num_atoms == None:
            for i in xrange(4):
                next(input_file)
            num_atoms = int(next(input_file).split()[0])
            while True:
                new_line = next(input_file)
                if new_line.find('--') == -1:
                    num_atoms = int(new_line.split()[0])
                else:
                    break
        elif (line.find("Atom  AN") > -1 or line.find("Atom AN") > -1) and num_atoms:
            vfreq1 = []
            vfreq2 = []
            vfreq3 = []
            for atom in xrange(num_atoms):
                reading_line = input_file.next()
                broken = reading_line.split()
                vfreq1.extend([float(broken[2]), float(broken[3]), float(broken[4])])
                vfreq2.extend([float(broken[5]), float(broken[6]), float(broken[7])])
                vfreq3.extend([float(broken[8]), float(broken[9]), float(broken[10])])
            modes.append(vfreq1)
            modes.append(vfreq2)
            modes.append(vfreq3)

    input_file.close()

    return modes


def read_frequencies(path):
    """
    Read normal modes

    path: string

    frequencies: list
    """
    input = open(path, "r")
    frequencies = []
    for line in input:
        if line.find(" Frequencies --") > -1:
            broken = line.split()
            f1 = float(broken[2])
            frequencies.append(f1)
            f2 = float(broken[3])
            frequencies.append(f2)
            f3 = float(broken[4])
            frequencies.append(f3)
    input.close()

    return frequencies


def read_molecule(path):
    pass


def create_prody_modes(frequencies, modes):
    my_modes = numpy.array(modes).T
    my_frequencies = numpy.array(frequencies)
    prody_modes = prody.NMA()
    prody_modes.setEigens(vectors=my_modes,values=my_frequencies)
    return prody_modes


def check_molecule(gaussian_molecule, gaudi_molecule):
    pass


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
        coords, elements, names, resnums, chids, betas, masses = [], [], [], [], [], [], []
        chimera2prody = {}
        offset_chimera_residue = min(r.id.position for r in molecule.residues)

        for i, atm in enumerate(molecule.atoms):
            chimera2prody[atm.coordIndex] = i
            coords.append(tuple(atm.coord()))  # array documentation to improve
            elements.append(atm.element.name)
            names.append(atm.name)
            resnums.append(atm.residue.id.position - offset_chimera_residue)
            chids.append(atm.residue.id.chainId)
            masses.append(atm.element.mass)
            betas.append(atm.bfactor)

        prody_molecule.setCoords(coords)
        prody_molecule.setElements(elements)
        prody_molecule.setNames(names)
        prody_molecule.setResnums(resnums)
        prody_molecule.setChids(chids)
        prody_molecule.setBetas(betas)
        prody_molecule.setMasses(masses)
        prody_molecule.setTitle(str(molecule.name))
        prody_molecule.setBonds([(chimera2prody[bond.atoms[0].coordIndex],
                                  chimera2prody[bond.atoms[1].coordIndex]) for bond in molecule.bonds])

    except AttributeError:
        raise TypeError('Attribute not found. Molecule must be a chimera.Molecule')

    return prody_molecule, chimera2prody


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
