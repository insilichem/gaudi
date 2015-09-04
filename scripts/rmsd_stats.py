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

import csv
import sys
import itertools
import numpy

with open(sys.argv[1]) as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    filtered_csv = filter(lambda x: float(x[2]) > 0, reader)
    grouped_csv = itertools.groupby(filtered_csv, lambda x: x[0])
    total, greats, goods, accs, avg_rmsd_first = 0, [], [], [], []
    for k, g in grouped_csv:
        group = list(g)
        pdb, result, rmsd = min(group, key=lambda x: x[2])
        rmsd = float(rmsd)
        avg_rmsd_first.append(next(
            float(c) for (a, b, c) in group if b == 'dock_000_Ligand.mol2'))
        total += 1
        if rmsd < 1.5:
            print(pdb, result, rmsd, 'GREAT!!')
            if result == 'dock_000_Ligand.mol2':
                accs.append(rmsd)
            greats.append(rmsd)
            goods.append(rmsd)
        elif rmsd < 2.5:
            print(pdb, result, rmsd, 'GOOD')
            goods.append(rmsd)
        else:
            print(pdb, result, rmsd)
    print('Poses < 1.5:', 100 * len(greats) / total,
          '%; Avg RMSD:', numpy.mean(greats))
    print('Poses < 2.5:', 100 * len(goods) / total,
          '%; Avg RMSD:', numpy.mean(goods))
    print('Accurate hits:', 100 * len(accs) / total,
          '%; Avg RMSD:', numpy.mean(accs))
    print('Avg RMSD of first result:', numpy.mean(avg_rmsd_first), 'A')
