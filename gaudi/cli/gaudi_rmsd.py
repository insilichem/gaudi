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
`gaudi.cli.gaudi_rmsd` calculates RMSD between two molecule poses in batch mode.
It should be provided with the path to the folders and the name of the
reference ligand. Run it with ``gaudi rmsd``.

.. note ::
    It requires rdkit for fast RMSD calculations.

"""

from rdkit import Chem
from rdkit.Chem import rdMolAlign
import argparse
import csv
import itertools
import numpy
import os
import yaml
import zipfile


def rmsd(dataset, outputfile, reference='reference.mol2', results='results'):
    """
    Take a benchmark essay and compute RMSD for every solution against the reference.
    The results are written to ``rmsd.csv`` in `dataset`.

    Parameters
    ----------
    dataset : str
        Path to previously benchmarked dataset
    outputfile : str
        The default name of the output file of each GAUDI essay
    reference : str, optional
        Name of the ligand reference molecule to calculate RMSD against
    results : str, optional
        Name of the results subdirectory

    Notes
    -----

    .. code-block:: console

        For each folder in dataset, open its results and:
            Extract ligand.mol2 from zipfile
            Open reference.mol2
            Calculate RMSD
        Print resulting RMSD values to CSV file
        If some error occurred, the RMSD value will be set to a negative float:
            -1.0: The reference molecule could not be loaded
            -2.0: The ligand molecule could not be loaded
            -3.0: Ref and ligand were loaded, but RMSD calculation failed
    """

    with open(os.path.join(dataset, 'rmsd.csv'), 'w+') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=['test', 'ligand', 'rmsd'])
        writer.writeheader()
        for d in os.listdir(dataset):
            try:
                reference = Chem.MolFromMol2File(os.path.join(dataset, d, reference),
                                                 sanitize=True, removeHs=False)
            except OSError:
                print("There was an error with reference from {}".format(d))
                writer.writerow(
                    {'test': d, 'ligand': 'REF_ERROR', 'rmsd': -1.0})
            else:
                for name, ligandblock in ligands_from_zip(dataset, d, results, outputfile):
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


def stats(dataset, great_threshold=1.5, good_threshold=2.5):
    """ 
    Run through the CSV file and get some stats and print outcome to 
    ``rmsd.txt`` in `dataset`.

    Parameters
    ----------
    dataset : str
        Path to previously benchmarked and rmsd'd dataset
    great_threshold : float
        If RMSD < this value, it's considered *great*
    great_threshold : float
        If RMSD < this value, it's considered *good*

    """

    with open(os.path.join(dataset, 'rmsd.csv')) as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        filtered_csv = filter(lambda x: float(x[2]) > 0, reader)
        grouped_csv = itertools.groupby(filtered_csv, lambda x: x[0])

    with open(os.path.join(dataset, 'rmsd.txt')) as txt:
        total, greats, goods, accs, avg_rmsd_first = 0, [], [], [], []
        for k, g in grouped_csv:
            group = list(g)
            pdb, result, rmsd = min(group, key=lambda x: x[2])
            rmsd = float(rmsd)
            avg_rmsd_first.append(next(
                float(c) for (a, b, c) in group if b.endswith('_000_Ligand.mol2')))
            total += 1
            if rmsd < great_threshold:
                txt.write('{}\t{}\t{}\tGREAT\n'.format(pdb, result, rmsd))
                if result.endswith('_000_Ligand.mol2'):
                    accs.append(rmsd)
                greats.append(rmsd)
                goods.append(rmsd)
            elif rmsd < good_threshold:
                txt.write('{}\t{}\t{}\GOOD\n'.format(pdb, result, rmsd))
                goods.append(rmsd)
            else:
                txt.write('{}\t{}\t{}\n'.format(pdb, result, rmsd))

        txt.write('Poses < 1.5: {}%\n\tAvg RMSD: {}A\n'.format(
                  100 * len(greats) / total, numpy.mean(greats)))
        txt.write('Poses < 2.5: {}%\n\tAvg RMSD: {}A\n'.format(
                  100 * len(goods) / total, numpy.mean(goods)))
        txt.write('Accurate hits: {}%\n\tAvg RMSD: {}A\n'.format(
                  100 * len(accs) / total, numpy.mean(accs)))
        txt.write('Avg RMSD of first result: {}A\n'.format(
                  numpy.mean(avg_rmsd_first)))


# Helpers
def ligands_from_zip(dataset, folder, results, outputfile):
    """
    Open all the listed solutions in ``*.gaudi-output`` file and create
    a generator that yields the name and mol2 data of the ligand of
    each solution.

    As a temporary workaround, the ligand is considered
    to be the lightest file with a mol2 extension.

    Parameters
    ----------
    dataset : str
        Path to benchmarked dataset being analyzed
    folder : str
        Path to the entry in dataset being analyzed
    results : str
        Path to the results subdirectory in each `folder`
    outputfile : str
        Name of the GAUDI output file in `results`

    """

    with open(os.path.join(dataset, folder, results, outputfile)) as f:
        data = yaml.load(f)
        for filename in data['GAUDI.results']:
            # Extract ligand.mol2 from zipfile
            # Open reference.mol2
            # Calculate RMSD
            # Print to file
            try:
                z = zipfile.ZipFile(
                    os.path.join(dataset, folder, results, filename))
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
    """
    Use RDKit molecule aligner just once, so it calculates the RMSD of both

    Parameters
    ----------
    ligand, reference : rdkit.Chem.molecule
    """

    if ligand is not None and reference is not None:
        return rdMolAlign.AlignMol(ligand, reference, maxIters=0)
    return -3.0


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Benchmarks a set given a template.')
    parser.add_argument('dataset', metavar='/path/to/dataset', type=str,
                        help='Location of the benchmark dataset')
    parser.add_argument('--outputfile', metavar='gaudi_file', type=str,
                        default='gaudi.gaudi-output', help='Name of GAUDI output file')
    parser.add_argument('--reference', metavar='ligand_name', type=str,
                        default='reference.mol2', help='Reference mol2 to benchmark against')
    parser.add_argument('--results', metavar='results_dir', type=str,
                        default='results', help='Subdirectory containing results')
    parser.add_argument('--good', type=float,
                        default=2.5, help='Subdirectory containing results')
    parser.add_argument('--great', type=float,
                        default=1.5, help='Subdirectory containing results')
    args = parser.parse_args()

    rmsd(args.dataset, args.inputfile, args.reference, args.results)
    stats(args.dataset, args.great, args.good)
