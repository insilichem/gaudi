#!/usr/bin/python

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
:mod:`gaudi.genes.molecule` implements a wrapper around Chimera.molecule objects
to expand its original features, such as appending new molecules.

This allows to build new structures with a couple of building blocks as a starting
point, as well as keeping several ligands as different potential solutions to the
essay (think about multi-molecule alternative docking). The user can also request
more :mod:`gaudi.genes.molecule` instances for the genome of the individual,
resulting in a competitive multi-docking essay.

To handle all this diversity, each construction is cached the first time is built.

This class is a dependency of most of the other genes (and even objectives), so it
will be requested almost always.

.. todo::

    As a result of the recent migration, some covalent interaction parameters
    are still present in the code. Ie, `vertex` argument in the keys. A solution
    must be proposed to get rid of those, probably a new gene that specifies this
    restrictions (`covalent`, for example).
"""

# Python
import os
import itertools
import random
import sys
import logging
# Chimera
import chimera
import BuildStructure
from chimera import UserError
from chimera.molEdit import addAtom, addBond
from WriteMol2 import writeMol2
# External dependencies
import deap
import yaml
from repoze.lru import LRUCache
# GAUDI
from gaudi import box, parse
from gaudi.genes import GeneProvider
from gaudi.genes import search

ZERO = chimera.Point(0.0, 0.0, 0.0)
logger = logging.getLogger(__name__)


def enable(**kwargs):
    return Molecule(**kwargs)


class Molecule(GeneProvider):

    """
    Interface around the :class:`gaudi.genes.molecule.Compound` to handle
    the GAUDI protocol and caching features.
    """
    _CATALOG = {}

    def __init__(self, path=None, symmetry=None,
                 **kwargs):
        GeneProvider.__init__(self, **kwargs)
        self._kwargs = kwargs
        self.path = path
        self.symmetry = symmetry
        self.origin = ZERO
        if 'gaudi.genes.search' in sys.modules:
            for g in self.parent.cfg.genes:
                if g.type == 'gaudi.genes.search' and g.target == self.name:
                    self.origin = g.origin
                    break

        try:
            self._compoundcache = self._cache['compounds']
        except KeyError:
            self._compoundcache = self._cache['compounds'] = LRUCache(300)
            self.catalog = self._CATALOG[
                self.name] = tuple(self._compile_catalog())

        self.catalog = self._CATALOG[self.name]
        self.allele = random.choice(self.catalog)
        self.compound = self.get(self.allele)

    def __deepcopy__(self, memo):
        """
        Fixes copy of Chimera unmutable objects (Atoms, Residues, etc)

        After each generation, DEAP will copy the current population using Python's
        __deepcopy___. Since Chimera uses some C++ wrappers for its own objects, such
        as atoms, residues and molecules, which don't implement __deepcopy___, this returns and error.
        These objects shouldn't be copied over anyway (that's why we have a cache), so overriding
        __deepcopy__ in this class circumvents the problem. All we do is ignore the recursivity
        of the original __deepcopy__ to prevent the access to the inner objects. All we keep is a
        reference to them.
        """
        new = self.__class__(self.path, self.symmetry,
                             **self._kwargs)
        new.__dict__.update((k, v) for k, v in self.__dict__.items())
        new.allele = self.allele + ()  # make sure we get a NEW allele
        return new

    def express(self):
        self.compound = self.get(self.allele)
        chimera.openModels.add([self.compound.mol], shareXform=True)
        box.pseudobond_to_bond(self.compound.mol)

        try:
            self.compound.place(chimera.Point(*self.origin))
        except TypeError:
            self.origin = search.parse_origin(self.origin, self.parent.genes)
            self.compound.place(chimera.Point(*self.origin))

    def unexpress(self):
        chimera.openModels.remove([self.compound.mol])
        # del self.compound

    def mate(self, mate):
        try:
            self.allele, mate.allele = deap.tools.cxTwoPoint(
                list(self.allele), list(mate.allele))
        except (StopIteration, ValueError):
            self.allele, mate.allele = mate.allele, self.allele
        else:
            self.allele, mate.allele = tuple(self.allele), tuple(mate.allele)

    def mutate(self, indpb):
        """
        VERY primitive. It only gets another compound.
        """
        if random.random() < self.indpb:
            self.allele = random.choice(self.catalog)

    def write(self, path, name):
        fullname = os.path.join(
            path, '{}_{}({}).mol2'.format(name, self.name, self.__class__.__name__))
        writeMol2([self.compound.mol], fullname, skip=[self.get_dummy_atoms()],
                  temporary=True, multimodelHandling='combined')
        return fullname

    ############
    def __getitem__(self, key):
        return self.get(key, 0)

    def get(self, key, vertex=0):
        # repoze.lru does not raise exceptions, so we must switch to LBYL
        compound = self._compoundcache.get((key, vertex))
        if not compound:
            compound = self.build(key)
            self._compoundcache.put((key, compound.vertex), compound)
        return compound

    def build(self, key, where=None):
        vertex = 0
        # if self.covalent:
        #   base = Compound()
        #   base.donor = base.add_dummy_atom(where.neighbors[0], serial=1)
        #   base.acceptor = base.add_dummy_atom(where, bonded_to=base.donor, serial=2)
        #   base.append(Compound(molecule=key[0], seed=random.random()))
        #   vertex = base.vertex
        # else:
        #   base = Compound(molecule=key[0])
        #   base.place(where)
        base = Compound(molecule=key[0])
        # base.place(self.origin)
        for molpath in key[1:]:
            base.append(Compound(molecule=molpath))

        # if self.flexible:
        #   base.update_rotatable_bonds()
        base.vertex = vertex
        return base

    def _compile_catalog(self):
        container = set()
        if os.path.isdir(self.path):
            folders = sorted(os.path.join(self.path, d) for d in os.listdir(self.path)
                             if os.path.isdir(os.path.join(self.path, d)) and not d.startswith('.')
                             and not d.startswith('_'))
            if folders:
                catalog = itertools.product(*[box.files_in(f, ext='mol2')
                                              for f in folders])
                if isinstance(self.symmetry, list):
                    folders_last_level = [os.path.basename(os.path.normpath(f))
                                          for f in folders]
                    for entry in catalog:
                        if all(os.path.basename(entry[folders_last_level.index(s1)]) ==
                                os.path.basename(
                                    entry[folders_last_level.index(s2)])
                                for (s1, s2) in self.symmetry):
                            self.catalog.append(entry)
                else:
                    container.update(tuple(catalog))
            else:
                container.update((f,)
                                 for f in box.files_in(self.path, ext='mol2'))
        elif os.path.isfile(self.path) and self.path.endswith('.mol2'):
            container.add((self.path,))
        return container

    def get_dummy_atoms(self, name='DUM'):
        for a in self.compound.mol.atoms:
            if a.name == name:
                yield a


class Compound(object):

    """
    Wraps `chimera.Molecule` instances and allows to perform copies,
    appending new fragments, free placement and extended attributes.

    .. todo::

        This was built a while a go (my first class), so it will
        probably need some refactoring.
    """

    def __init__(self, molecule=None, origin=None, seed=0.0, **kwargs):
        if isinstance(molecule, chimera.Molecule):
            self.mol = molecule
        elif not molecule or molecule == 'dummy':
            self.mol = _dummy_mol('dummy')
        else:
            self.mol, = chimera.openModels.open(molecule)
            chimera.openModels.remove([self.mol])

        self.mol.gaudi = self
        self.rotatable_bonds = []
        self.built_atoms = []
        self.donor = self.mol.atoms[0]
        self.acceptor = self.mol.atoms[-1]
        self.parse_attr()
        self.origin = origin
        self.seed = seed
        self.vertex = 0
        for k, v in kwargs.items():
            self.__dict__[k] = v
        if not hasattr(self, 'nonrotatable'):
            self.nonrotatable = []

    def parse_attr(self):
        try:
            f = open(self.mol.openedAs[0][:-4] + 'attr')
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

    def update_attr(self, d):
        for k, v in self.__dict__.items():
            if k in ('mol', 'angles', 'built_atoms'):
                continue
            if isinstance(v, list):
                setattr(self, k, [d[v_] if v_ in d else v_ for v_ in v])
            else:
                setattr(self, k, d[v] if v in d else v)

    def destroy(self):
        chimera.openModels.close([self.mol])
        del self

    ####
    # Building and moving functions
    ####
    def add_dummy_atom(self, where, name='dum', element=None, residue=None,
                       bonded_to=None, serial=None):
        if isinstance(where, chimera.Atom):
            element = where.element if not element else element
            where = where.coord()
        else:
            element = chimera.Element('C')

        residue = self.mol.residues[-1] if not residue else residue
        return addAtom(name, element, residue, where, serial, bonded_to)

    def append(self, molecule):
        self.attach(molecule, self.acceptor, molecule.donor)

    def prepend(self, molecule):
        self.attach(molecule, self.donor, molecule.donor)

    def attach(self, molecule, acceptor, donor):
        if acceptor not in self.mol.atoms:
            raise UserError('Specified atom is not part of molecule.')

        molecule.place_for_bonding(acceptor)
        self.vertex = molecule.vertex
        if not donor:
            donor = molecule.donor
        updated_atoms = self.join(molecule, acceptor, donor)

        # Update Gaudi ATTR
        self.update_attr(updated_atoms)  # update existant atoms
        self.acceptor = updated_atoms[molecule.acceptor]  # change acceptor
        self.axis_end = updated_atoms[molecule.axis_end]
        if hasattr(molecule, 'nonrotatable'):
            nonrot_atoms = [updated_atoms[a] for a in molecule.nonrotatable]
            if hasattr(self, 'nonrotatable'):
                self.nonrotatable.extend(nonrot_atoms)
            else:
                self.nonrotatable = nonrot_atoms

        molecule.destroy()

    def join(self, molecule, acceptor, donor, newres=False):
        target = acceptor
        sprouts = [donor]
        res = _dummy_mol().residues[0] if newres else target.residue
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

        self.built_atoms.append({a.serialNumber: b for (a, b) in built_atoms.items()})
        return built_atoms

    def place(self, where, anchor=None):
        if isinstance(where, chimera.Atom):
            where = where.coord()
        if not anchor:
            anchor = self.donor
        search.translate(self.mol, anchor, where)

    def place_for_bonding(self, target, anchor=None, seed=None):
        if not isinstance(target, chimera.Atom):
            raise UserError('Target must be a chimera.Atom object.')
        if not anchor:
            anchor = self.donor
        if not seed:
            seed = self.seed
        # Get target position
        target_pos, self.vertex = _new_atom_position(
            target, anchor.element, seed)
        # Place it
        self.place(target_pos)
        # Fix orientation
        anchor_pos, i = _new_atom_position(anchor, target.element)
        search.rotate(
            self.mol, [target.coord(), anchor.coord(), anchor_pos], 0.0)


def _dummy_mol(name):
    m = chimera.Molecule()
    m.name = name
    r = m.newResidue(name, 'het', 1, ' ')
    r.isHet = True
    return m


def _new_atom_position(atom, newelement, seed=0.0):
    try:
        geometry = chimera.idatm.typeInfo[atom.idatmType].geometry
    except KeyError:
        geometry = 3
        logger.warning("Using %s geometry for atom %s in molecule %s",
                       geometry, atom, molecule.name)
    bond_length = chimera.Element.bondLength(atom.element, newelement)
    neighbors_crd = [a.coord() for a in atom.neighbors]
    points = chimera.bondGeom.bondPositions(atom.coord(), geometry, bond_length,
                                            neighbors_crd)
    return points[int(seed * len(points))], int(seed * len(points))
