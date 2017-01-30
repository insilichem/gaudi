#!/usr/bin/env python
# -*- coding: utf-8 -*-

##############
# GaudiMM: Genetic Algorithms with Unrestricted
# Descriptors for Intuitive Molecular Modeling
# Authors:  Jaime Rodriguez-Guerra Pedregal
#            <jaime.rodriguezguerra@uab.cat>
#           Jean-Didier Marechal
#            <jeandidier.marechal@uab.cat>
# Web: https://bitbucket.org/insilichem/gaudi
##############

"""
`gaudi.cli.gaudi_benchmark` enables sequential runs of a given GAUDI essay 
in all folders contained under given path. Useful for benchmarks. 
Run it with ``gaudi benchmark``.

"""

from __future__ import print_function
import os
import yaml
from gaudi.cli import gaudi_run

def main(dataset, templatefile):
    """
    Helper script to benchmark a dataset with a base GAUDI input.

    Load a common template with at least two ``molecule`` instances:
    one called ``Protein`` and a second called ``Ligand``. Then, it replaces
    the original paths in the template with those found in the subsequent
    directories and launches a new Chimera instance to perform the essay.

    Parameters
    ----------
    dataset : str
        Path to the benchmark dataset
    templatefile : str
        Path to the base GAUDI input
    """

    # Load input template
    with open(templatefile, 'r') as f:
        template = yaml.load(f)

    # Find indices of Protein and Ligand genes
    prot_index = next(i for (i, g) in enumerate(template['genes']) if g['name'] == 'Protein')
    ligand_index = next(i for (i, g) in enumerate(template['genes']) if g['name'] == 'Ligand')

    protein = template['genes'][prot_index]['path']
    ligand = template['genes'][ligand_index]['path']
    results = template['output']['path']

    try:
        os.mkdir(results)
    except (OSError, IOError):
        pass

    # Scan directories in supplied dataset
    dataset = os.path.abspath(dataset)
    dirs = os.listdir(dataset)
    for i, d in enumerate(dirs):
        # 1 - Perform calculation
        print("Processing item {}: {}".format(i+1, d))
        # Replace original output path with a subfolder in current dataset
        template['output']['path'] = os.path.join(results, d)
        template['output']['name'] = d

        # Replace original Ligand and Protein paths with new ones
        template['genes'][prot_index]['path'] = os.path.join(dataset, d, protein)
        template['genes'][ligand_index]['path'] = os.path.join(dataset, d, ligand)

        # Write the overwritten template to a temp location
        input_filename = os.path.join(results, d + '.yaml')
        with open(input_filename, 'w') as f:
            f.write(yaml.dump(template))

        # Launch the essay!
        gaudi_run.main(input_filename, debug=False)

        # Report progress after job is completed
        print("Processed {}/{} ({}%)".format(i+1, len(dirs), 100 * (i+1) / len(dirs)))
