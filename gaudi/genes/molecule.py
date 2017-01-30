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
This gene implements a wrapper around Chimera.molecule objects
to expand its original features, such as appending new molecules.

This allows to build new structures with a couple of building blocks as a starting
point, as well as keeping several ligands as different potential solutions to the
essay (think about multi-molecule alternative docking). The user can also request
more :mod:`gaudi.genes.molecule` instances for the genome of the individual,
resulting in a competitive multi-docking essay.

To handle all this diversity, each construction is cached the first time is built.

This class is a dependency of most of the other genes (and even objectives), so it
will be requested almost always.

"""

# Python
from __future__ import print_function
import os
import itertools
import random
import logging
import tempfile
try:
    from cStringIO import StringIO
except ImportError:
    from StringIO import StringIO
# Chimera
import chimera
from _multiscale import get_atom_coordinates
from chimera import UserError
from chimera.molEdit import addAtom, addBond
from AddH import simpleAddHydrogens
from WriteMol2 import writeMol2
# External dependencies
import deap
import yaml
from boltons.cacheutils import LRU
from pdbfixer import PDBFixer
from simtk.openmm.app import PDBFile
# GAUDI
from gaudi import box, parse
from gaudi.genes import GeneProvider
from gaudi.genes import search
from gaudi.exceptions import AtomsNotFound, ResiduesNotFound, TooManyAtoms, TooManyResidues


ZERO = chimera.Point(0.0, 0.0, 0.0)
logger = logging.getLogger(__name__)


def enable(**kwargs):
    kwargs = Molecule.validate(kwargs)
    return Molecule(**kwargs)


class Molecule(GeneProvider):

    """
    Interface around the :class:`gaudi.genes.molecule.Compound` to handle
    the GAUDI protocol and caching features.

    Parameters
    ----------
    path : str, optional
        Path to a molecule file or a directory containing dirs of molecule files.

    symmetry : str, optional
        If `path` is a directory, list of pairs of directories whose chosen
        mocule must be the same, thus enabling *symmetry*.

    hydrogens : bool, optional
        Add hydrogens to Molecule (True) or not (False).

    pdbfix : bool, optional
        Only for testing and debugging. Better run pdbfixer prior to GAUDI.
        Fix potential issues that may cause troubles with OpenMM forcefields.

    Attributes
    ----------
    _CATALOG : dict
        Class attribute (shared among all `Molecule` instances) that holds
        all the possible molecules GAUDI can build given current `path`. 

        If `path` is a single molecule file, that's the only possibility, but
        if it's set to a directory, the engine can potentially build all the
        combinations of molecule blocks found in those subdirectories.

        Normally, it is accessed via the `catalog` property.

    Notes
    -----
    **Use of `_cache`**

    `Molecule` class uses `_cache` to store already built molecules. Its entry in
    the `_cache` directory is a `boltons.cacheutils.LRU` cache (least recently used)
    set to a maximum size of 300 entries.

    .. todo ::

        The LRU cache size should be proportional to essay size, depending on
        the population size and number of generations, but also taking available
        memory into account (?).


    """

    _validate = {
        parse.Required('path'): parse.RelPathToInputFile(),
        'symmetry': [[basestring]],
        'hydrogens': parse.Boolean,
        'pdbfix': parse.Boolean,
        }

    _CATALOG = {}
    SUPPORTED_FILETYPES = ('mol2', 'pdb')

    def __init__(self, path=None, symmetry=None, hydrogens=False, pdbfix=False, **kwargs):
        self._kwargs = kwargs.copy()
        GeneProvider.__init__(self, **kwargs)
        self._kwargs = kwargs
        self.path = path
        self.symmetry = symmetry
        self.hydrogens = hydrogens
        self.pdbfix = pdbfix
        try:
            self.catalog = self._CATALOG[self.name]
        except KeyError:
            self.catalog = self._CATALOG[self.name] = tuple(self._compile_catalog())
        self._compounds_cache = self._cache.setdefault(self.name + '_compounds', LRU(300))
        self._atomlookup_cache = self._cache.setdefault(self.name + '_atomlookup', LRU(300))
        self._residuelookup_cache = self._cache.setdefault(self.name + '_residuelookup', LRU(300))
        self.allele = random.choice(self.catalog)

        # An optimization for similarity methods: xform coords are
        # cached here after all genes have expressed. See Individual.express.
        self._expressed_coordinates = None

    def __deepcopy__(self, memo):
        new = self.__class__(path=self.path, symmetry=self.symmetry, hydrogens=self.hydrogens,
                             pdbfix=self.pdbfix, **self._kwargs)
        new.allele = self.allele + ()
        new._expressed_coordinates = self._expressed_coordinates
        return new

    @property
    def compound(self):
        """
        Get expressed allele on-demand (read-only attribute)
        """
        return self.get(self.allele)

    def express(self):
        """
        Adds Chimera molecule object to the viewer canvas.

        It also converts pseudobonds (used by Chimera to depict
        coordinated ligands to metals) to regular bonds.
        """
        chimera.openModels.add([self.compound.mol], shareXform=True)
        box.pseudobond_to_bond(self.compound.mol)

    def unexpress(self):
        """
        Removes the Chimera molecule from the viewer canvas
        (without deleting the object).
        """
        chimera.openModels.remove([self.compound.mol])

    def mate(self, mate):
        """
        .. todo::

            Allow mating while preserving symmetry
        """
        if len(self.allele) == 1:
            return
        if not self.symmetry:
            try:
                self.allele, mate.allele = deap.tools.cxTwoPoint(
                    list(self.allele), list(mate.allele))
            except (StopIteration, ValueError):
                self.allele, mate.allele = mate.allele, self.allele
            else:
                self.allele, mate.allele = tuple(
                    self.allele), tuple(mate.allele)

    def mutate(self, indpb):
        """
        VERY primitive. It only gets another compound.
        """
        if random.random() < self.indpb:
            self.allele = random.choice(self.catalog)

    def write(self, path=None, name=None, absolute=None, combined_with=None, filetype='mol2'):
        """
        Writes full mol2 to disk.

        .. todo::

            It'd be preferable to get a string instead of a file
        """
        if path and name:
            fullname = os.path.join(
                path, '{}_{}.{}'.format(name, self.name, filetype))
        elif absolute:
            fullname = absolute
        else:
            fileobject, fullname = tempfile.mkstemp(prefix='gaudi', suffix='.{}'.format(filetype))
            logger.warning("No output path provided. Using tempfile %s.", fullname)

        molecules = [self.compound.mol]
        if combined_with:
            molecules.extend(ind.compound.mol for ind in combined_with)
        if filetype == 'mol2':
            writeMol2(molecules, fullname, temporary=True, multimodelHandling='combined')
        elif filetype == 'pdb':
            chimera.pdbWrite(molecules, self.compound.mol.openState.xform, fullname)
        else:
            raise ValueError('Filetype {} not recognized. Try with mol2 or pdb'.format(filetype))
        return fullname

    @classmethod
    def clear_cache(cls):
        GeneProvider.clear_cache()
        cls._CATALOG.clear()

    ############
    def __getitem__(self, key):
        """Implements dict-like item retrieval"""
        return self.get(key)

    def get(self, key):
        """
        Looks for the compound corresponding to `key` in `_cache`. If found,
        return it. Else, build it on demand, store it on cache and return it.

        Parameters
        ----------
        key : str
            Path (or combination of) to the requested molecule. It should be
            extracted from `catalog`.

        Returns
        -------
        gaudi.genes.molecule.Compound
            The result of building the requested molecule.

        """
        try:
            return self._compounds_cache[key]
        except KeyError:
            self._compounds_cache[key] = compound = self.build(key)
            return compound

    def build(self, key, where=None):
        """
        Builds a `Compound` following the recipe contained in `key` through
        `Compound.append` methods.

        Parameters
        ----------
        key : tuple of str
            Paths to the molecule blocks that comprise the final molecule.
            A single molecule is just a one-block recipe.

        Returns
        -------
        Compound
            The final molecule result of the sequential appending.
        """
        base = Compound(molecule=key[0])
        for molpath in key[1:]:
            base.append(Compound(molecule=molpath))
        if self.hydrogens:
            base.add_hydrogens()
        if self.pdbfix:
            base.apply_pdbfix()
        return base

    def _compile_catalog(self):
        """
        Computes all the possible combinations of given directories and mol2 files,
        taking symmetry requirements into account.

        Parameters
        ----------
        self.path
        self.symmetry

        Returns
        -------
        set
            A set of tuples of paths, each indicating a build recipe for a molecule
        """
        container = set()
        if os.path.isdir(self.path):
            folders = sorted(os.path.join(self.path, d) for d in os.listdir(self.path)
                             if os.path.isdir(os.path.join(self.path, d)) and not d.startswith('.')
                             and not d.startswith('_'))
            if folders:
                catalog = itertools.product(*[box.files_in(f, ext=self.SUPPORTED_FILETYPES)
                                              for f in folders])
                if isinstance(self.symmetry, list):
                    folders_last_level = [os.path.basename(os.path.normpath(f))
                                          for f in folders]
                    for entry in catalog:
                        if all(os.path.basename(entry[folders_last_level.index(s1)]) ==
                                os.path.basename(entry[folders_last_level.index(s2)])
                                for (s1, s2) in self.symmetry):
                            container.add(entry)
                else:
                    container.update(tuple(catalog))
            else:
                container.update((f,) for f in box.files_in(self.path,
                                                            ext=self.SUPPORTED_FILETYPES))
        elif (self.path.split('.')[-1] in self.SUPPORTED_FILETYPES):
            container.add((self.path,))
        return container

    # API methods
    def find_atoms(self, serial, only_one=False):
        if serial == '*':
            return self.compound.mol.atoms
        try:
            return self._atomlookup_cache[(self.allele, serial)]
        except KeyError:
            atoms = self._find_atoms(serial, only_one)
            self._atomlookup_cache[(self.allele, serial)] = atoms
            return atoms

    def _find_atoms(self, serial, only_one=False):
        if isinstance(serial, int):  # search by serial number
            atoms = [a for a in self.compound.mol.atoms if a.serialNumber == serial]
        elif isinstance(serial, basestring):  # search by name
            try:
                atoms = [getattr(self.compound, serial)]
            except AttributeError:
                atoms = [a for a in self.compound.mol.atoms if a.name == serial]
        else:
            raise ValueError('Serial {} not valid'.format(serial))
        if atoms:
            if only_one and len(atoms) > 1:
                raise TooManyAtoms("Found {} atoms for serial {} but expected 1.".format(len(atoms), serial))
            return atoms
        raise AtomsNotFound("Atom '{}' not found in {}".format(serial, self.name))

    def find_atom(self, serial):
        return self.find_atoms(serial, only_one=True)[0]

    def find_residues(self, position, only_one=False):
        if position == '*':
            return self.compound.mol.residues
        try:
            return self._residuelookup_cache[(self.allele, position)]
        except KeyError:
            residues = self._find_residues(position, only_one)
            self._residuelookup_cache[(self.allele, position)] = residues
            return residues

    def _find_residues(self, position, only_one):
        residues = [r for r in self.compound.mol.residues if r.id.position == position]
        if residues:
            if only_one and len(residues) > 1:
                raise TooManyResidues("Found {} residues for position {}")
            return residues
        raise ResiduesNotFound("Residue {} not found in {}".format(position, self.name))

    def find_residue(self, position):
        return self.find_residues(position, only_one=True)[0]

    def xyz(self, transformed=True):
        return get_atom_coordinates(self.compound.mol.atoms, transformed=transformed)


class Compound(object):

    """
    Wraps `chimera.Molecule` instances and allows to perform copies,
    appending new fragments, free placement and extended attributes.

    Parameters
    ----------
    molecule : chimera.Molecule or str, optional
        The chimera.Molecule object to wrap, or

        'dummy' (for a empty molecule), or

        The path to a mol2 file that will be parsed by Chimera to return 
        a valid chimera.Molecule object.

    origin : 3-tuple of float, optional
        Coordinates to the place where the molecule should be placed

    seed : float, optional
        Random seed, used by the vertex chooser on appending.

    kwargs, optional
        Additional parameters that should be injected as attributes

    Attributes
    ----------
    mol : chimera.Molecule
    donor : chimera.Atom
        If self were to be bonded to another Compound, this atom would be
        the one involved in the bond. By default, first element in
        chimera.Molecule.atoms
    acceptor : chimera.Atom
        If another Compound wanted to be bonded to self, this atom would
        be the one involved in the bond. By default, last element in
        chimera.Molecule.atoms
    origin : 3-tuple of float
        Default location of `mol` in the 3D canvas.
    seed : float
        Randomness seed for vertex computation.
    rotatable_bonds : list of chimera.Bond
        Bonds that can be torsioned
    nonrotatable : list of chimera.Bond
        Bonds that, although possible, should not be torsiond. Ie, fixed bonds.
    built_atoms : list of chimera.Atom
        Memo of already built_atoms

    Notes
    -----
    .. todo::

        Instead of a `molecule` parameter overload, think of using
        equivalent classmethods.

    """

    def __init__(self, molecule=None, origin=None, seed=0.0, **kwargs):
        if isinstance(molecule, chimera.Molecule):
            self.mol = molecule
        elif not molecule or molecule == 'dummy':
            self.mol = _dummy_mol('dummy')
        else:
            molecules = chimera.openModels.open(molecule)
            self.mol = molecules[0]
            chimera.openModels.remove(molecules)

        self.mol.gaudi = self
        self.rotatable_bonds = []
        self.built_atoms = []
        self.donor = self.mol.atoms[0]
        self.acceptor = self.mol.atoms[-1]
        self.parse_attr()
        self.origin = origin
        self.seed = seed
        for k, v in kwargs.items():
            self.__dict__[k] = v
        if not hasattr(self, 'nonrotatable'):
            self.nonrotatable = []

    def parse_attr(self):
        try:
            basepath, _ = os.path.splitext(self.mol.openedAs[0])
            f = open(basepath + '.attr')
        except IOError:
            logger.warning("No attr file found for molecule %s",
                           self.mol.openedAs[0])
        else:
            attr = yaml.load(f)
            flat_attr = {}
            if 'atoms' in attr:
                for k, v in attr['atoms'].items():
                    flat_attr[k] = next(
                        a for a in self.mol.atoms if a.serialNumber == v)
            if 'angles' in attr:
                flat_attr['angles'] = {}
                for k, v in attr['angles'].items():
                    flat_attr['angles'][tuple(k)] = v
            if 'bonds' in attr:
                for k, v in attr['bonds'].items():
                    if k == 'nonrotatable':
                        if isinstance(v, list):
                            flat_attr[k] = box.atoms_by_serial(
                                *v, atoms=self.mol.atoms)
                        elif v.strip().lower() in ('all', 'yes'):
                            flat_attr[k] = self.mol.atoms

            self.__dict__.update(flat_attr)
            f.close()

    def update_attr(self, d):
        for k, v in self.__dict__.items():
            if k in ('mol', 'angles', 'built_atoms'):
                continue
            if isinstance(v, list):
                setattr(self, k, [d[v_] if v_ in d else v_ for v_ in v])
            else:
                setattr(self, k, d[v] if v in d else v)

    def destroy(self):
        """
        Removes `mol` from canvas, deletes it and then, deletes itself
        """
        chimera.openModels.close([self.mol])
        del self

    ####
    # Building and moving functions
    ####
    def add_dummy_atom(self, where, name='dum', element=None, residue=None,
                       bonded_to=None, serial=None):
        """
        Adds a placeholder atom at the coordinates specified by `where`

        Parameters
        ----------
        where : chimera.Atom or 3-tuple of float
            Coordinates of target location. A chimera.Atom can be supplied,
            in which case its coordinates will be used (via `.coord()`)
        name : str, optional
            Name for the new atom
        element : chimera.Element, optional
            Element of the new atom
        residue : chimera.Residue, optional
            Residue that will incorporate the new atom
        bonded_to : chimera.Atom, optional
            Atom that will form a bond with new atom
        serial : int
            Serial number that will be assigned to atom

        """
        if isinstance(where, chimera.Atom):
            element = where.element if not element else element
            where = where.coord()
        else:
            element = chimera.Element('C')

        residue = self.mol.residues[-1] if not residue else residue
        return addAtom(name, element, residue, where, serial, bonded_to)

    def append(self, molecule):
        """
        Wrapper around `attach` to add a new molecule to `self`, 
        using `self.acceptor` as bonding atom.

        Parameters
        ----------
        molecule : Compound
        """
        self.attach(molecule, self.acceptor, molecule.donor)

    def prepend(self, molecule):
        """
        Wrapper around `attach` to add a new molecule to `self`, 
        using `self.donor` as bonding atom.

        Parameters
        ----------
        molecule : Compound
        """
        self.attach(molecule, self.donor, molecule.donor)

    def attach(self, molecule, acceptor, donor):
        """
        Call `join` to bond `molecule` to `self.mol` and updates attributes.

        Parameters
        ----------
        molecule : Compound
            The molecule that will be attached to `self`
        acceptor : chimera.Atom
            Atom of `self` that will participate in the new bond
        donor : chimera.Atom
            Atom of `molecule` that will participate in the new bond

        Notes
        -----
        .. note ::

            After joining the two molecules together, we have to update the 
            attributes of self to match the new molecular reality. For example,
            the atom participating in the bond will not be available for new 
            bonds (we are forcing linear joining for now), so the new `donor`
            will be inherited from `molecule`, before deleting the object.

        """
        if acceptor not in self.mol.atoms:
            raise UserError('Specified atom is not part of molecule.')

        molecule.place_for_bonding(acceptor)
        if not donor:
            donor = molecule.donor
        updated_atoms = self.join(molecule, acceptor, donor)

        # Update Gaudi ATTR
        self.update_attr(updated_atoms)  # update existant atoms
        self.acceptor = updated_atoms[molecule.acceptor]  # change acceptor
        self.axis_end = updated_atoms[molecule.axis_end]
        if hasattr(molecule, 'flagged'):
            self.flagged = updated_atoms[molecule.flagged]
        if hasattr(molecule, 'nonrotatable'):
            nonrot_atoms = [updated_atoms[a] for a in molecule.nonrotatable]
            if hasattr(self, 'nonrotatable'):
                self.nonrotatable.extend(nonrot_atoms)
            else:
                self.nonrotatable = nonrot_atoms

        molecule.destroy()

    def join(self, molecule, acceptor, donor, newres=False):
        """
        Take a molecule and bond it to `self.mol`, copying the atoms in the 
        process so they're contained in the same `chimera.Molecule` object.

        Parameters
        ----------
        molecule : Compound
            The molecule that will be attached to `self`
        acceptor : chimera.Atom
            Atom of `self` that will participate in the new bond
        donor : chimera.Atom
            Atom of `molecule` that will participate in the new bond
        newres : bool
            If True, don't reuse `acceptor.residue` for the new molecule,
            and create a new blank one instead.

        Returns
        -------
        built_atoms : dict
            Maps original atoms in `molecules` to their new counterparts
            in `self.mol`.

        Notes
        -----
        .. note ::

            Chimera does not allow bonds between different chimera.Molecule
            objects, so firstly, we have to copy the atoms of ``molecule`` to 
            ``self.mol`` and, only then, make the joining bond.

            It traverses the atoms of ``molecule`` and adds a copy of each
            of them to ``self.mol`` using ``chimera.molEdit.addAtom`` in the same spot
            of space. All the bonds are preserved and, finally, bond the two molecules.

            The algorithm starts by adding the bonding atom of ``molecule`` (``donor``), to
            the ``sprouts`` list. Then, the loop starts:

            .. code-block:: python

                while sprouts contains atoms:
                    sprout = sprouts.pop(0)
                    copy sprout to self.mol
                    for each neighbor of sprout
                        copy neighbor to self.mol
                        if neighbor itself has more than one neighbor (ie, sprout)
                            add neighbor to sprouts 

            Along the way, we have to take care of already processed sprouts, so we
            don't repeat ourselves. That's what the built_atoms dict is for.

            Also, instead of letting addAtom guess new serial numbers, we calculate
            them beforehand by computing the highest serial number in ``self.mol`` 
            prior to the additions and then incremeting one by one on a per-element
            basis.

        .. todo ::

            This code is UGLY. Find a better way!
        """
        target = acceptor
        sprouts = [donor]
        res = _dummy_mol().residues[0] if newres else target.residue
        # I think I was trying to keep a copy of the original residue,
        # but this does not copy anything. Since lists are mutable, we
        # are just aliasing res.atoms with res_atoms.
        # Maybe I want something like:
        # res_atoms = res.atoms[:] # ?
        res_atoms = res.atoms

        i = max(a.serialNumber for a in self.mol.atoms)
        index = box.highest_atom_indices(self.mol)
        for a in molecule.mol.atoms:
            try:
                index[a.element.name] += 1
            except KeyError:
                index[a.element.name] = 1
            finally:
                a.name = a.element.name + str(index[a.element.name])
        built_atoms = {}
        while sprouts:
            sprout = sprouts.pop(0)  # get first atom
            if sprout.name in res.atomsMap:
                target = res.atomsMap[sprout.name][-1]
            else:
                i += 1
                built_atoms[sprout] = target = addAtom(sprout.name,
                                                       sprout.element, res, sprout.coord(),
                                                       bondedTo=target, serialNumber=i)

            for a in sprout.neighbors:
                if a.element.number == 1:
                    continue
                if a.name not in res.atomsMap:
                    needBuild = True
                else:
                    # atom is already present, but it can be part of a cycle
                    # if we get to it it's because another atom is linking it
                    needBuild = False
                    built = res.atomsMap[a.name][-1]
                if needBuild:
                    i += 1
                    built_atoms[a] = built = addAtom(a.name, a.element, res, a.coord(),
                                                     bondedTo=target, serialNumber=i)
                    # if a has more than one neighbor:
                    if len(a.neighbors) > 1:
                        sprouts.append(a)  # this new atom can be a new sprout
                # link!
                if built not in target.bondsMap and built not in res_atoms:
                    addBond(target, built)

        self.built_atoms.append(
            {a.serialNumber: b for (a, b) in built_atoms.items()})
        return built_atoms

    def place(self, where, anchor=None):
        """
        Convenience wrapper around `translate`, supplying default

        Parameters
        ----------
        where : 3-tuple of float, or chimera.Atom
            Coordinates of destination. If chimera.Atom, use supplied `coord()`

        anchor : chimera.Atom
            The atom that will guide the translation.

        Notes
        -----
        .. note ::
            *Anchor atoms* are called that way in the API because I picture them
            as the one we pick with our hands to drag the molecule to the desired
            place.
        """
        if isinstance(where, chimera.Atom):
            where = where.coord()
        if not anchor:
            anchor = self.donor
        search.translate(self.mol, anchor, where)

    def place_for_bonding(self, target, anchor=None):
        """
        Translate `self.mol` to a colavent distance of `target` atom,
        with an adequate orientation.

        Parameters
        ----------
        target : chimera.Atom
            The atom we are willing to bond later on.

        anchor : chimera.Atom, optional

        """
        if not isinstance(target, chimera.Atom):
            raise UserError('Target must be a chimera.Atom object.')
        if not anchor:
            anchor = self.donor
        # Get target position
        target_pos = _new_atom_position(target, anchor.element)
        # Place it
        self.place(target_pos)
        # Fix orientation
        anchor_pos = _new_atom_position(anchor, target.element)
        search.rotate(self.mol, [target.coord(), anchor.coord(), anchor_pos], 0.0)

    def add_hydrogens(self):
        """
        Add missing hydrogens to current molecule
        """
        simpleAddHydrogens([self.mol])

    def apply_pdbfix(self, pH=7.0):
        """
        Run PDBFixer and replace original molecule with new one
        """
        self.mol = _apply_pdbfix(self.mol, pH)


def _apply_pdbfix(molecule, pH=7.0, add_hydrogens=False):
    """
    Run PDBFixer to ammend potential issues in PDB format.

    Parameters
    ----------
    molecule : chimera.Molecule
        Chimera Molecule object to fix.
    pH : float, optional
        Target pH for adding missing hydrogens.
    add_hydrogens : bool, optional
        Whether to add missing hydrogens or not.

    Returns
    -------
    memfile : StringIO
        An in-memory file with the modified PDB contents
    """
    memfile = StringIO()
    chimera.pdbWrite([molecule], chimera.Xform(), memfile)
    chimera.openModels.close([molecule])
    memfile.seek(0)
    fixer = PDBFixer(pdbfile=memfile)
    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.removeHeterogens(True)
    if add_hydrogens:
        fixer.addMissingHydrogens(pH)
    memfile.close()

    memfile = StringIO()
    PDBFile.writeFile(fixer.topology, fixer.positions, memfile)
    memfile.seek(0)
    molecule = chimera.openModels.open(memfile, type="PDB", identifyAs=molecule.name)
    chimera.openModels.remove(molecule)
    memfile.close()
    return molecule[0]


def _dummy_mol(name='dummy'):
    """
    Create an empty molecule, with an empty residue. Used for new molecules
    we intend to create on a per-atom basis.

    Parameters
    ----------
    name : str, optional

    Returns
    -------
    chimera.Molecule

    """
    m = chimera.Molecule()
    m.name = name
    r = m.newResidue(name, 'het', 1, ' ')
    r.isHet = True
    return m


def _new_atom_position(atom, newelement, seed=0.0):
    """
    Get suitable coordinates for new atoms that intend to bond to `atom`.

    Parameters
    ----------
    atom : chimera.Atom
        Atom new atoms will be trying to bond to
    newelement : chimera.Element
        Element of the atom that will be bonding `atom`.
    seed : float, optional
        Residual kwarg of old API. REMOVE THIS!

    Returns
    -------
    chimera.Point
        The first coordinate set in the returned list of possible points.

    Notes
    -----
    .. note ::
        `newelement` is needed because the bond length depends on the
        atom types involved in such bond. Ie, C-C != C-N.
    """
    try:
        geometry = chimera.idatm.typeInfo[atom.idatmType].geometry
    except KeyError:
        geometry = 3
        logger.warning("Using %s geometry for atom %s in molecule %s",
                       geometry, atom, atom.molecule.name)
    bond_length = chimera.Element.bondLength(atom.element, newelement)
    neighbors_crd = [a.coord() for a in atom.neighbors]
    points = chimera.bondGeom.bondPositions(atom.coord(), geometry, bond_length,
                                            neighbors_crd)
    return points[0]
