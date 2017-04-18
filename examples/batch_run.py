#!/usr/bin/env pychimera
# -*- coding:utf-8 -*-

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

from __future__ import print_function, absolute_import, division
import os
import sys
from copy import deepcopy
from multiprocessing import Pool, cpu_count

import chimera
from WriteMol2 import writeMol2
from SplitMolecule.split import split_molecules
from gaudi.cli.gaudi_run import main as gaudi_run
from gaudi.parse import Settings


def stdout_to_file(workspace):
    path = os.path.join(workspace, 'log')
    path = _incremental_existing_path(path, separator="__")
    sys.stdout = open(path + '.stdout', "w")
    sys.stderr = open(path + '.stderr', "w")

def _incremental_existing_path(path, separator="__"):
    keep_trying = 1
    while os.path.exists(path):
        base, ext = os.path.splitext(path)
        base = base.rsplit(separator, 1)[0]
        path = '{}{}{}{}'.format(base, separator, keep_trying, ext)
        keep_trying += 1
    return path

def prepare_workspace(name, *molecules):
    basedir, filename = os.path.split(name)
    workspace_name, _ = os.path.splitext(filename)
    workspace_path = os.path.join(basedir, 'benchmark', workspace_name)
    workspace_path = _incremental_existing_path(workspace_path)
    os.makedirs(workspace_path)

    paths = []
    for i, molecule in enumerate(molecules):
        molname = '{}_{}.mol2'.format(os.path.splitext(molecule.name)[0].replace(' ', '_'), i)
        mol2file = os.path.join(workspace_path, molname)
        writeMol2([molecule], mol2file, temporary=True)
        paths.append(mol2file)

    return workspace_path, paths


def split_metal_protein(molecule):
    """
    Take a molecule file and separate the metal ions from the
    rest of the molecule to different files.

    Parameters
    ----------
    molecule : chimera.Molecule
        Molecule to edit

    Returns
    -------
    protein, metal : chimera.Molecule
    """
    # Retrieve proteins and write them to single mol2 file
    metals = [a for a in molecule.atoms if a.element in chimera.elements.metals]
    split_molecules(molecules=[molecule], atoms=[metals])
    return sorted(chimera.openModels.list(), key=lambda m: m.numAtoms)


def list_molecules(path):
    """
    List molecule-containing files in path
    """
    for name in os.listdir(path):
        _, ext = os.path.splitext(name)
        if ext in ('.mol2', '.pdb') and os.path.isfile(name):
            yield os.path.join(path, name)


def clean_canvas(molecule):
    chimera.runCommand('delete solvent')
    for element in ('CA', 'NA'):
        chimera.runCommand('sel @/element={}'.format(element))
        if len(chimera.selection.currentAtoms()) > 1:
            chimera.runCommand('del @/element={}'.format(element))
    chimera.runCommand('delete @/element=NA')


def benchmark(args):
    cfg, molfile = args
    molecule = chimera.openModels.open(molfile)[0]
    clean_canvas(molecule)
    splitted = split_metal_protein(molecule)
    if len(splitted) != 2:
        print("! File {} does not contain metals... Skipping!".format(molfile),
              file=sys.stderr)
        return
    metals, protein = splitted
    workspace, paths = prepare_workspace(molfile, protein, metals)
    prot_gene = next(g for g in cfg.genes if g.name == 'Protein')
    metal_gene = next(g for g in cfg.genes if g.name == 'Metal')
    prot_gene.path, metal_gene.path = paths[0:2]
    cfg.output.path = os.path.join(workspace, 'output')
    cfg.output.name = os.path.splitext(molfile)[0]
    chimera.openModels.close(chimera.openModels.list())
    stdout_to_file(workspace)
    cfg.validate()
    try:
        print("Running GAUDI for {}".format(cfg.output.name))
        gaudi_run(cfg)
    except KeyboardInterrupt:
        raise Exception('Ctrl+C!')

def main(template, path, n_processes=cpu_count()):
    if n_processes > cpu_count() or n_processes < 1:
        n_processes = cpu_count()
    pool = Pool(processes=n_processes, maxtasksperchild=1)
    cfg = Settings(template, validation=False)
    print("Running {} at {} with {} processes...".format(template, path, n_processes))
    try:
        args = [(deepcopy(cfg), molfile) for molfile in list_molecules(path)]
        pool.map(benchmark, args, chunksize=1)
    except KeyboardInterrupt:
        pool.terminate()
    except Exception as e:
        print(e)
        pool.terminate()
    finally:
        pool.close()
        pool.join()

if __name__ == '__main__':
    try:
        main(sys.argv[1], sys.argv[2], *map(int, sys.argv[3:4]))
    except IndexError:
        sys.exit('Usage: pychimera batch_run.py '
                 '/path/to/template.yaml /path/to/benchmark/files [number of processes]')