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
This script calculates RMSD between two molecule poses in batch mode.
It should be provided with the path to the folders and the name of the
reference ligand.
"""

import argparse
import csv
import os
from rdkit import Chem
from rdkit.Chem import rdMolAlign
import yaml
import zipfile


def parse_args():
    parser = argparse.ArgumentParser(
        description='Benchmarks a set given a template.')
    parser.add_argument('dataset', metavar='/path/to/dataset', type=str,
                        help='Location of the benchmark dataset.')
    parser.add_argument('--file', metavar='gaudi_file', type=str,
                        default='gaudi.yaml', help='Name of GAUDI output file')
    parser.add_argument('--ligand', metavar='ligand_name', type=str,
                                            default='ligand.mol2', help='Tagname of the ligand file')
    parser.add_argument('--results', metavar='results_dir', type=str,
                        default='results', help='Subdirectory containing results')

    return parser.parse_args()


def ligands_from_zip(d, args):
    with open(os.path.join(args.dataset, d, args.results, args.file)) as f:
        data = yaml.load(f)
        for filename in data['GAUDI.results']:
            # Extract ligand.mol2 from zipfile
            # Open reference.mol2
            # Calculate RMSD
            # Print to file
            try:
                z = zipfile.ZipFile(
                    os.path.join(args.dataset, d, args.results, filename))
            except zipfile.BadZipfile:
                print("Not a valid GAUDI result")
            else:
                mol2details = [info for info in z.infolist()
                               if info.filename.endswith('.mol2')]
                mol2details.sort(key=lambda i: i.file_size)
                yield mol2details[0].filename, z.open(mol2details[0].filename).read()
            finally:
                z.close()


def calculate_rmsd(ligand, reference):
    if ligand is not None and reference is not None:
        return rdMolAlign.AlignMol(ligand, reference, maxIters=0)
    return -3.0


def main():
    args = parse_args()
    with open('rmsd2.csv', 'w+') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=['test', 'ligand', 'rmsd'])
        writer.writeheader()
        for d in os.listdir(args.dataset):
            ref_path = os.path.join(args.dataset, d, args.ligand)
            try:
                reference = Chem.MolFromMol2File(ref_path,
                                                 sanitize=True, removeHs=False)
            except OSError:
                print("There was an error with reference from {}".format(d))
                writer.writerow(
                    {'test': d, 'ligand': 'REF_ERROR', 'rmsd': -1.0})
            else:
                for name, ligandblock in ligands_from_zip(d, args):
                    try:
                        ligand = Chem.MolFromMol2Block(ligandblock,
                                                       sanitize=True, removeHs=False)
                    except OSError:
                        print(
                            "There was an error with ligand {} from {}".format(name, d))
                        writer.writerow(
                            {'test': d, 'ligand': name, 'rmsd': -2.0})
                    else:
                        rmsd = calculate_rmsd(ligand, reference)
                        writer.writerow(
                            {'test': d, 'ligand': name, 'rmsd': rmsd})

if __name__ == '__main__':
    main()
