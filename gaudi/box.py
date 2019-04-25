#!/usr/bin/env python
# -*- coding: utf-8 -*-

##############
# GaudiMM: Genetic Algorithms with Unrestricted
# Descriptors for Intuitive Molecular Modeling
#
# https://github.com/insilichem/gaudi
#
# Copyright 2017 Jaime Rodriguez-Guerra, Jean-Didier Marechal
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
##############

"""
This module is a messy collection of useful functions used all along GAUDI.

.. todo::

    Some of these functions are hardly used, so maybe we should clean it a
    little in the future...

"""

# Python
import os
import cProfile
from contextlib import contextmanager
import sys
# Chimera
import chimera
import numpy as np
from scipy.sparse.csgraph import shortest_path as _shortest_path
# Gaudi
import gaudi


def atoms_between(atom1, atom2):
    """
    Finds all connected atoms between two given atoms. Can be slow with
    large molecules.
    """
    chain1 = [atom1]
    chain2 = [atom2]
    i = 0
    while i < len(chain1):
        a1 = chain1[i]
        if atom2 not in a1.neighbors:
            chain1.extend([a for a in a1.neighbors if a not in chain1])
        i += 1
    j = 0
    while j < len(chain2):
        a2 = chain2[j]
        if atom1 not in a2.neighbors:
            chain2.extend([a for a in a2.neighbors if a not in chain2])
        j += 1
    return set(chain1) & set(chain2)


def atoms_by_serial(*serials, **kw):
    """
    Find atoms in kw['atoms'] with serialNumber = `serials`.

    Parameters
    ----------
    serials : int
        List of serial numbers to match
    atoms : list of chimera.Atom, optional
        List of atoms to be traversed while looking for serial numbers

    Returns
    -------
    list of chimera.Atom

    """
    if not kw['atoms']:
        kw['atoms'] = [a for m in chimera.openModels.list() for a in m.atoms]
    return [a for a in kw['atoms'] if a.serialNumber in serials]


def create_single_individual(path):
    """
    Create an individual within Chimera. Convenience method for Chimera IDLE.
    """
    def prepare_input(path, parser):
        """
        Parses input file and validate paths
        """
        import os
        import logging

        def build_path(basedir, path):
            """
            Processes tildes and join paths to base directory of input file.
            ``os.path.join`` is smart enough to not join two absolute paths, returning
            the last one provided. ``os.path.normpath`` simplifies joined paths by
            parsing residual dots or double dots.
            """
            return os.path.normpath(os.path.join(basedir, os.path.expanduser(path)))

        # Parse input
        try:
            # os.path.realpath prepends the working directory to relative paths
            path = os.path.abspath(os.path.expanduser(path))
        except IndexError:
            print "ERROR: Input file not provided."
        else:
            cfg = parser(path)
            inputdir = os.path.dirname(path)

        # Tilde expansion in paths and abs/rel path support
        cfg.general.outputpath = build_path(inputdir, cfg.general.outputpath)
        for g in cfg.genes:
            if g.module == 'gaudi.genes.molecule':
                g.path = build_path(inputdir, g.path)
                if not os.path.exists(g.path):
                    print "ERROR: Path " + g.path + " is wrong. Check your input file.\n"

        # Create dirs
        try:
            os.makedirs(cfg.general.outputpath)
        except OSError:
            if os.path.isfile(cfg.general.outputpath):
                print "ERROR: Output path is already a file. Please change it.\n"

        # Register loggers and handlers for both stdout and file
        logger = logging.getLogger('gaudi')
        logger.setLevel(logging.DEBUG)

        # create CONSOLE handler and set level to error
        handler = logging.StreamHandler()
        handler.setLevel(logging.DEBUG)
        formatter = logging.Formatter("%(levelname)s - %(message)s")
        handler.setFormatter(formatter)
        logger.addHandler(handler)

        return cfg

    import deap
    from deap import creator, tools
    from . import base, plugin, parse
    cfg = prepare_input(path, parse.Settings)
    plugin.import_plugins(*cfg.genes)
    plugin.import_plugins(*cfg.objectives)
    toolbox = deap.base.Toolbox()
    toolbox.register("call", (lambda fn, *args, **kwargs: fn(*args, **kwargs)))
    toolbox.register("individual", toolbox.call, base.Individual, cfg)
    ind = toolbox.individual()
    environment = gaudi.base.Environment(cfg)
    return ind, environment


def do_cprofile(func):
    """
    Decorator to cProfile a certain function and output the results
    to cprofile.out
    """
    def profiled_func(*args, **kwargs):
        profile = cProfile.Profile()
        try:
            profile.enable()
            result = func(*args, **kwargs)
            profile.disable()
            return result
        finally:
            profile.dump_stats('cprofile.out')
    return profiled_func


def draw_interactions(interactions, startCol='FF0000', endCol='FFFF00',
                      key=None, name="Custom pseudobonds"):
    """
    Draw pseudobonds depicting atoms relationships.

    Parameters
    ----------
    interactions : list of tuples
        Each tuple contains an interaction, defined, at least,
        by the two atoms involved.
    startCol : str, optional
        Hex code for the initial color of the pseudobond
        (closer to the first atom of the pair).
    endCol : str, optional
        Hex code for the final color of the pseudobond.
        (closer to the second atom of the pair)
    key : int, optional
        The index of an interaction tuple that represent the alpha
        channel in the color used to depict the interaction.
    name : str, optional
        Name of the pseudobond group created.

    Returns
    -------
    chimera.pseudoBondGroup

    """
    if not len(interactions):
        return
    pb = chimera.misc.getPseudoBondGroup(name)
    color = _hex_to_rgb(startCol) + [1.0]
    if key:
        max_ = max(abs(_[3]) for _ in interactions)
    for i in interactions:
        npb = pb.newPseudoBond(i[0], i[1])
        if key is not None:
            intensity = (max_ - abs(i[key])) / (max_)
            opacity = 1 - 0.7 * intensity
            if startCol != endCol:
                color = _linear_color(intensity, startCol, endCol) + [opacity]
            else:
                color = _hex_to_rgb(startCol) + [opacity]
        npb.color = chimera.MaterialColor(*color)
    return pb


def files_in(path, ext=None):
    """
    Returns all the files in a given directory, filtered by extension if desired.

    Parameters
    ----------
    path : str
    ext : list of str, optional
        File extension(s) to filter on.

    Returns
    -------
        List of absolute paths
    """
    if ext:
        if isinstance(ext, str):
            ext = [ext]
        return [os.path.join(path, fn) for fn in next(os.walk(path))[2]
                if fn.split('.')[-1] in ext]
    return [os.path.join(path, fn) for fn in next(os.walk(path))[2]]


def find_nearest(anchor, atoms):
    """
    Find the atom of `atoms` that is closer to `anchor`, in terms of
    number of atoms in between.

    Note
    ----
    Can be very slow with lots of atoms in a molecule. For large
    structures, using ``gaudi.box.shortest_distance_matrix`` as
    done in ``gaudi.genes.torsion.Torsion._compute_rotatable_bonds``
    is recommended.
    """
    try:
        return next(a for a in atoms if a is anchor)
    except StopIteration:  # oops, didn't find it...
        return min(atoms, key=lambda a: len(atoms_between(anchor, a)))


def highest_atom_indices(r):
    """
    Returns a dictionary with highest atom indices in given residue
    Key: value -> element.name: highest index in residue
    """
    results = {'C': 0, 'H': 0, 'N': 0, 'O': 0}
    for a in r.atoms:
        if a.name[1:].isdigit():
            atom = a.name[:1]
            num = int(a.name[1:])
            if atom not in results:
                results[atom] = num
            elif results[atom] < num:
                results[atom] = num

        elif a.name[2:].isdigit():
            atom = a.name[:2]
            num = int(a.name[2:])
            if atom not in results:
                results[atom] = num
            elif results[atom] < num:
                results[atom] = num
    return results


def incremental_existing_path(path, separator="__"):
    keep_trying = 1
    while os.path.exists(path):
        base, ext = os.path.splitext(path)
        base = base.rsplit(separator, 1)[0]
        path = '{}{}{}{}'.format(base, separator, keep_trying, ext)
        keep_trying += 1
    return path
_incremental_existing_path = incremental_existing_path


@contextmanager
def open_models_and_close(*args, **kwargs):
    models = chimera.openModels.open(*args, **kwargs)
    yield models
    chimera.openModels.close(models)


def pseudobond_to_bond(molecule, remove=False):
    """
    Transforms every pseudobond in `molecule` to a covalent bond

    Parameters
    ----------
    molecule : chimera.Molecule
    remove : bool
        If True, remove original pseudobonds after actual bonds
        have been created.

    """
    pbgroup = chimera.misc.getPseudoBondGroup(
        "coordination complexes of %s (%s)" %
        (molecule.name, molecule))  # , associateWith=[molecule])
    if pbgroup.pseudoBonds:
        for pb in pbgroup.pseudoBonds:
            chimera.molEdit.addBond(*pb.atoms)
            if remove:
                pbgroup.deletePseudoBond(pb)
        pbm = molecule.pseudoBondMgr()
        pbm.deletePseudoBondGroup(pbgroup)


def shortest_distance_matrix(molecule, graph=None, atom_map=None):
    """
    Compute minimum distance matrix for all atoms in molecule.

    Parameters
    ----------
    molecule : chimera.Molecule

    graph : array_like
        Matrix-form graph with NxN elements (``N=molecule.numAtoms``)
        , where ``[i,j]`` values are 1 if atom ``i`` is bonded to ``j``.
        If not provided, it will be built.

    atom_map : dict, optional
        Mapping from ``chimera.Atom`` objects to its position
        in ``molecule.atoms``. If not provided, it will be built.

    Returns
    -------
    distances : array_like
        Specific distances can be queried in ``distances[i,j]``,
        with ``i,j`` being the indices of each atom in ``molecule.atoms``
    graph
    atom_map
    """
    if not atom_map:
        atom_map = {a: i for (i, a) in enumerate(molecule.atoms)}
    if not graph:
        graph = np.zeros((molecule.numAtoms, molecule.numAtoms))
        for b in molecule.bonds:
            i1, i2 = atom_map[b.atoms[0]], atom_map[b.atoms[1]]
            graph[i1, i2] = 1
    distances = _shortest_path(graph, method='auto', directed=False,
                               unweighted=True)
    return distances, graph, atom_map


def suppress_ksdssp(trig_name, my_data, molecules):
    """
    Monkey-patch Chimera triggers to disable KSDSSP computation
    """
    for m in molecules.created:
        m.structureAssigned = True


def write_individuals(inds, outpath, name, evalfn, remove=True):
    """
    Write an individual to disk.

    .. note ::

        Deprecated since an Individual object is able to write itself
        to disk.
    """
    from WriteMol2 import writeMol2
    if not os.path.isdir(outpath):
        os.makedirs(outpath)
    objectivesnames = (o for o in inds[0].fitness.objectives)
    header = ' '.join('{:>10}'.format(x) for x in objectivesnames)
    results = ['{:>{len_}} {}'.format('Filename', header, len_=len(name) + 10)]
    for i, ind in enumerate(inds):
        ind.express()
        ligand = ind.genes['Ligand'].compound
        fullname = '{}{}__{:03d}.mol2'.format(outpath, name, i + 1)
        writeMol2([ligand.mol], fullname, temporary=True,
                  multimodelHandling='combined')

        fitness = ' '.join('{:>10.6g}'.format(x) for x in ind.fitness.values)
        rotamers = ind.parsed_rotamers if 'rotamers' in ind.genes else []
        with open(fullname, 'r+') as f:
            mol2data = f.readlines()
            mol2data[:0] = ['#~ Generated by GAUDI\n\n']
            mol2data[-1:] = ['@<TRIPOS>COMMENT\n',
                             'GAUDI.score\n',
                             header, '\n',
                             fitness,
                             '\n/GAUDI.score\n',
                             '\nGAUDI.rotamers\n'] + \
                rotamers + \
                ['\n/GAUDI.rotamers']
            f.seek(0)
            f.write(''.join(mol2data))

        if remove:
            chimera.openModels.remove([ligand.mol])
        results.append('{} {}'.format(fullname.split('/')[-1], fitness))
    return results


def sequential_bonds(atoms, s):
    """
    Returns bonds in `atoms` in sequential order, beginning at atom `s`
    """
    if s not in atoms:
        atoms.append(s)
    bonds = list(
        set([b for a in atoms for b in a.bonds if set(b.atoms) <= set(atoms)]))
    nbonds = []
    while bonds:
        b = bonds.pop(0)
        if s in b.atoms:
            nbonds.append(b)
            s = b.otherAtom(s)
        else:
            bonds.append(b)
    return nbonds


def rmsd(a, b):
    import math
    if isinstance(a, chimera.Molecule):
        a = a.atoms
    if isinstance(b, chimera.Molecule):
        b = b.atoms

    a.sort(key=lambda z: z.serialNumber)
    b.sort(key=lambda z: z.serialNumber)
    sqdist = sum(x.xformCoord().sqdistance(y.xformCoord())
                 for x, y in zip(a, b))

    return math.sqrt(sqdist / float(len(a)))


@contextmanager
def silent_stdout():
    _stdout = sys.stdout
    with open(os.devnull, 'w') as sys.stdout:
        yield
    sys.stdout = _stdout


def stdout_to_file(workspace, stderr=True):
    path = os.path.join(workspace, 'log')
    path = _incremental_existing_path(path)
    sys.stdout = open(path + '.stdout', "w")
    if stderr:
        sys.stderr = open(path + '.stderr', "w")


def _hex_to_rgb(hexa):
    return [int(hexa[i:i + 2], 16) for i in range(0, 6, 2)]


def _linear_color(value, start, end):
    color = []
    c = 0
    for s, e in zip(_hex_to_rgb(start), _hex_to_rgb(end)):
        s, e = sorted([s, e])
        c = s + value * (e - s)
        color.append(c / 255.)
    return color
