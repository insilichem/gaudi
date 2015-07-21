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
Document this!
"""

import argparse
import os
import subprocess
import sys
import tempfile
import yaml

gaudi_launch = os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])),
                            'launch.py')

parser = argparse.ArgumentParser(
    description='Benchmarks a set given a template.')

parser.add_argument('template', metavar='/path/to/gaudi/template', type=str,
                    help='A GAUDI input file that will serve as template.')
parser.add_argument('dataset', metavar='/path/to/dataset', type=str,
                    help='Location of the benchmark dataset.')
parser.add_argument('--protein', metavar='protein_name', type=str,
                    default='protein.mol2', help='Tagname of the protein file')
parser.add_argument('--ligand', metavar='ligand_name', type=str,
                                        default='ligand.mol2', help='Tagname of the ligand file')
parser.add_argument('--chimera', metavar='/path/to/bin/chimera', type=str,
                    default='~/.local/UCSF-Chimera64-1.10.1/bin/chimera',
                    help='Path to UCSF Chimera binaries')
parser.add_argument('--gaudi', metavar='/path/to/gaudi/launch.py', type=str,
                    default=gaudi_launch)
args = parser.parse_args()

with open(args.template, 'r') as f:
    template = yaml.load(f)

dirs = os.listdir(args.dataset)
for i, d in enumerate(dirs):
    print "Processing item", d
    template['general']['outputpath'] = os.path.join(args.dataset, d,
                                                     'results')
    template['genes'][0]['path'] = os.path.join(args.dataset, d, args.protein)
    template['genes'][1]['path'] = os.path.join(args.dataset, d, args.ligand)

    number, filename = tempfile.mkstemp()
    with open(filename, 'w') as f:
        f.write(yaml.dump(template))

    command = [os.path.expanduser(args.chimera), '--nogui',
               '--silent', '--script', '{} {}'.format(args.gaudi, filename)]
    subprocess.call(command)
    print "Processed {}/{} ({}%)".format(i, len(dirs), 100 * i / len(dirs))
