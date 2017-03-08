#!/usr/bin/env pychimera
# -*- coding:utf-8 -*-

##############
# GaudiMM: Genetic Algorithms with Unrestricted
# Descriptors for Intuitive Molecular Modeling
# 
# http://bitbucket.org/insilichem/gaudi
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
from multiprocessing import Pool

import chimera
from WriteMol2 import writeMol2
from SplitMolecule.split import split_molecules
from gaudi.cli.gaudi_run import main as gaudi_run
from gaudi.parse import Settings

def stdout_to_file():
    sys.stdout = open(str(os.getpid()) + ".out", "w")

def prepare_workspace(name, *molecules):
    basedir, filename = os.path.split(name)
    workspace_name, _ = os.path.splitext(filename)
    workspace_path = os.path.join(basedir, 'benchmark', workspace_name)
    
    keep_trying = 1
    while keep_trying:
        try:
            os.makedirs(workspace_path)
        except OSError:
            if os.path.exists(workspace_path):
                workspace_path = '{}___{}'.format(workspace_path.rsplit('___', 1)[0], keep_trying)
                keep_trying += 1
                continue
        else:
            keep_trying = 0

    paths = []
    for i, molecule in enumerate(molecules):
        molname = '{}_{}.mol2'.format(os.path.splitext(molecule.name)[0].replace(' ', '_'), i)
        mol2file = os.path.join(workspace_path, molname)
        writeMol2([molecule], mol2file, temporary=True)
        paths.append(mol2file)

    return workspace_path, paths

def split_metal_protein(molecule):
    """
    Take a molecule file and extract the metals and protein
    in separate files.

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
    return sorted(chimera.openModels.list(), key= lambda m: m.numAtoms)


def list_molecules(path):
    """
    List molecule-containing files in path
    """
    for name in os.listdir(path):
        _, ext = os.path.splitext(name)
        if ext in ('.mol2', '.pdb') and os.path.isfile(name):
            yield os.path.join(path, name)

def clean_canvas():
    chimera.runCommand('delete solvent')
    chimera.runCommand('delete @/element=CA')
    chimera.runCommand('delete @/element=CL')
    chimera.runCommand('delete @/element=NA')

def benchmark(args):
    cfg, molfile = args
    molecule = chimera.openModels.open(molfile)[0]
    clean_canvas()
    metals, protein = split_metal_protein(molecule)
    workspace, paths = prepare_workspace(molfile, protein, metals)
    prot_gene = next(g for g in cfg.genes if g.name == 'Protein')
    metal_gene = next(g for g in cfg.genes if g.name == 'Metal')
    prot_gene.path = paths[0]
    metal_gene.path = paths[1]
    cfg.output.path = os.path.join(workspace, 'output')
    cfg.output.name = os.path.splitext(molfile)[0]
    chimera.openModels.close(chimera.openModels.list())
    cfg.validate(cfg)
    try:
        gaudi_run(cfg)
    except KeyboardInterrupt:
        raise Exception('Ctrl+C!')

def main(template, path):
    pool = Pool(initializer=stdout_to_file)
    cfg = Settings(template, validation=False)
    try:
        pool.map(benchmark, [(deepcopy(cfg), molfile) for molfile in list_molecules(path)])
    except KeyboardInterrupt:
        pool.terminate()
    except Exception as e:
        print(e)
        pool.terminate()


if __name__ == '__main__':
    try:
        main(sys.argv[1], sys.argv[2])
    except IndexError:
        sys.exit('Usage: batch_run.py /path/to/template.yaml /path/to/benchmark/files')