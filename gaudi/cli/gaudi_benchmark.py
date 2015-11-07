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
Sequential runs of a given GAUDI essay in all folders contained under
given path. Useful for benchmarks.

"""

from __future__ import print_function
import argparse
import os
import subprocess
import sys
import tempfile
import yaml


def arguments():
    gaudi_launch = os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])),
                                'launch.py')
    parser = argparse.ArgumentParser(
        description='Benchmarks a set given a template.')
    parser.add_argument('template', metavar='/path/to/template.gaudi-input', type=str,
                        help='A GAUDI input file that will serve as template.')
    parser.add_argument('dataset', metavar='/path/to/dataset', type=str,
                        help='Location of the benchmark dataset.')
    parser.add_argument('--protein', metavar='protein_name', type=str,
                        default='protein.mol2', help='Filename of the protein file')
    parser.add_argument('--ligand', metavar='ligand_name', type=str,
                        default='ligand.mol2', help='Filename of the ligand file. '
                        'It is supposed to be the real solution, as well.')
    parser.add_argument('--results', metavar='folder_name', type=str,
                        default='results', help='Name of the output folder.')
    parser.add_argument('--chimera', metavar='/path/to/bin/chimera', type=str,
                        default='~/.local/UCSF-Chimera64-1.10.1/bin/chimera',
                        help='Path to UCSF Chimera binary')
    parser.add_argument('--gaudi', metavar='/path/to/gaudi/launch.py', type=str,
                        default=gaudi_launch)

    return parser.parse_args()


def run(args):
    """
    Load a common template with at least two ``molecule`` instances:
    one called ``Protein`` and a second called ``Ligand``. Then, it replaces
    the original paths in the template with those found in the subsequent 
    directories and launches a new Chimera instance to perform the essay.
    """

    # Load input template
    with open(args.template, 'r') as f:
        template = yaml.load(f)

    # Find indices of Protein and Ligand genes
    prot_index = next(
        i for (i, g) in enumerate(template['genes']) if g['name'] == 'Protein')
    ligand_index = next(
        i for (i, g) in enumerate(template['genes']) if g['name'] == 'Ligand')

    # Scan directories in supplied dataset
    dirs = os.listdir(args.dataset)
    for i, d in enumerate(dirs):
        # 1 - Perform calculation
        print("Processing item ", i+1, ": ", d)
        # Replace original output path with a subfolder in current dataset
        template['general']['outputpath'] = os.path.join(args.dataset, d,
                                                         args.results)

        # Replace original Ligand and Protein paths with new ones
        template['genes'][prot_index]['path'] = os.path.join(
            args.dataset, d, args.protein)
        template['genes'][ligand_index]['path'] = os.path.join(
            args.dataset, d, args.ligand)

        # Write the overwritten template to a temp location
        number, filename = tempfile.mkstemp()
        with open(filename, 'w') as f:
            f.write(yaml.dump(template))

        # Launch the essay!
        command = [os.path.expanduser(args.chimera), '--nogui',
                   '--silent', '--script', '{} {}'.format(args.gaudi, filename)]
        subprocess.call(command)

        # Report progress after job is completed
        print("Processed {}/{} ({}%)".format(
            i+1, len(dirs), 100 * (i+1) / len(dirs)))


def main():
    args = arguments()
    run(args)

if __name__ == "__main__":
    main()
