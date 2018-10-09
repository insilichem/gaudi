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
This objective is a wrapper around the scoring functions provided by
`CCDC's GOLD <https://www.ccdc.cam.ac.uk/solutions/csd-discovery/components/gold/>`_.

It will use the rescoring abilities in GOLD to extract the fitness corresponding to
any of the available scoring functions:

- GoldScore (goldscore)
- ChemScore (chemscore)
- Astex Statistical Potential (asp)
- CHEMPLP (chemplp)

Since GOLD is commercial software, you will need to install it separately
and provide a valid license! This is just a wrapper. Make sure to set
all the needed environment variables, such as CCDC_LICENSE_FILE, and
that 'gold_auto' is in $PATH. Check tests/test_objectives_gold.py for an
example; make sure to have GOLDXX/bin before GOLDXX/GOLD/bin!

"""

# Python
from distutils.spawn import find_executable
from tempfile import _get_default_tempdir as default_tempdir, _get_candidate_names as tempnames
from string import Template
import logging
import numpy as np
import os
import subprocess
import sys
# Chimera
from Molecule import atom_positions
from WriteMol2 import writeMol2
# GAUDI
from gaudi.objectives import ObjectiveProvider
from gaudi import parse

logger = logging.getLogger(__name__)


def enable(**kwargs):
    kwargs = Gold.validate(kwargs)
    return Gold(**kwargs)


_TEMPLATE = Template("""
  GOLD CONFIGURATION FILE

  AUTOMATIC SETTINGS
autoscale = 1

  POPULATION
popsiz = auto
select_pressure = auto
n_islands = auto
maxops = auto
niche_siz = auto

  GENETIC OPERATORS
pt_crosswt = auto
allele_mutatewt = auto
migratewt = auto

  FLOOD FILL
radius = $RADIUS
origin = $ORIGIN
do_cavity = 0
floodfill_atom_no = 0
cavity_file =
floodfill_center = point

  DATA FILES
ligand_data_file $LIGAND 1
param_file = DEFAULT
set_ligand_atom_types = 1
set_protein_atom_types = 0
directory = .
tordist_file = DEFAULT
make_subdirs = 0
save_lone_pairs = 0
fit_points_file = fit_pts.mol2
read_fitpts = 0

  FLAGS
internal_ligand_h_bonds = 0
flip_free_corners = 0
match_ring_templates = 0
flip_amide_bonds = 0
flip_planar_n = 0
flip_pyramidal_n = 0
rotate_carboxylic_oh = fix
use_tordist = 1
postprocess_bonds = 0
fix_all_protein_rotatable_bonds = 1
solvate_all = 1

  TERMINATION
early_termination = 1
n_top_solutions = 3
rms_tolerance = 1.5

  CONSTRAINTS
force_constraints = 0

  COVALENT BONDING
covalent = 0

  SAVE OPTIONS
save_score_in_file = 1
save_protein_torsions = 0
clean_up_option delete_all_initialised_ligands
clean_up_option delete_rank_file

  FITNESS FUNCTION SETTINGS
initial_virtual_pt_match_max = 3
relative_ligand_energy = 0
gold_fitfunc_path = $SCORING
score_param_file = DEFAULT

  RUN TYPE
run_flag = RESCORE no_simplex no_file no_strip

  PROTEIN DATA
protein_datafile = $PROTEIN

""")


class Gold(ObjectiveProvider):
    """
    Gold class

    Parameters
    ----------
    protein : str
        The name of molecule acting as protein
    ligand : str
        The name of molecule acting as ligand
    scoring : str, optional, defaults to chemscore
        Fitness function to use. Choose between chemscore, chemplp,
        goldscore and asp.
    score_component : str, optional, defaults to 'Score'
        Scoring fields to parse out of the rescore.log file, such as
        Score, DG, S(metal), etc.
    radius : float, optional, defaults to 10.0
        Radius (in A) of binding site sphere, the origin of which is
        automatically centered at the ligand's center of mass.

    Returns
    -------
    float
        Interaction energy as reported by GOLD's chosen scoring function
    """
    _validate = {
        parse.Required('protein'): parse.Molecule_name,
        parse.Required('ligand'): parse.Molecule_name,
        'scoring': parse.In(['chemscore', 'chemplp', 'goldscore', 'asp']),
        'radius': parse.Coerce(float),
        'score_component': str,
    }

    def __init__(self, protein='Protein', ligand='Ligand',
                 scoring='chemscore', score_component='Score',
                 radius=10, *args, **kwargs):
        ObjectiveProvider.__init__(self, **kwargs)
        self.protein_names = [protein]
        self.ligand_names = [ligand]
        self.scoring = scoring
        self.score_component = score_component
        self.radius = radius
        self.executable = find_executable('gold_auto')
        if self.executable is None:
            sys.exit('GOLD could not be found in $PATH. Is it (correctly) installed?')
        self._oldworkingdir = os.getcwd()
        self._paths = {}
        if os.name == 'posix' and os.path.exists('/dev/shm'):
            self.tmpdir = '/dev/shm'
        else:
            self.tmpdir = default_tempdir()

    def get_molecule_by_name(self, ind, *names):
        """
        Get a molecule gene instance of individual by its name
        """
        for name in names:
            yield ind.find_molecule(name)

    def evaluate(self, ind):
        """
        Run a subprocess calling LigScore binary with provided options,
        and parse the results. Clean tmp files at exit.
        """
        self.tmpfile = os.path.join(self.tmpdir, next(tempnames()))
        proteins = list(self.get_molecule_by_name(ind, *self.protein_names))
        ligands = list(self.get_molecule_by_name(ind, *self.ligand_names))

        protein_path = self.prepare_proteins(proteins)
        ligand_path = self.prepare_ligands(ligands)
        origin = self.origin(ligands[0])
        command = self.prepare_command(protein_path, ligand_path, origin)

        try:
            os.chdir(self.tmpdir)
            p = subprocess.call(command)
            return self.parse_output('rescore.log')
        except (subprocess.CalledProcessError, IOError):
            logger.warning("Could not run GOLD with command %s", command)
            return -100000 * self.weight
        finally:
            self.clean()
            os.chdir(self._oldworkingdir)

    def prepare_proteins(self, proteins):
        proteinpath = '{}_proteins.pdb'.format(self.tmpfile)
        last_protein = proteins.pop()
        last_protein.write(
            absolute=proteinpath, combined_with=proteins, filetype='pdb')
        self._paths['proteins'] = proteinpath
        return proteinpath

    def prepare_ligands(self, ligands):
        ligandpath = '{}_ligands.mol2'.format(self.tmpfile)
        ligand_mols = [lig.compound.mol for lig in ligands]

        writeMol2(ligand_mols, ligandpath, temporary=True, multimodelHandling='combined')
        self._paths['ligands'] = ligandpath
        return ligandpath

    def origin(self, molecule):
        molecule = molecule.compound.mol
        coordinates = atom_positions(molecule.atoms, molecule.openState.xform)
        masses = np.fromiter((a.element.mass for a in molecule.atoms),
                             dtype='float32', count=molecule.numAtoms)
        return np.average(coordinates, axis=0, weights=masses)

    def prepare_command(self, protein_path, ligand_path, origin):
        replaces = dict(PROTEIN=protein_path, LIGAND=ligand_path,
                        ORIGIN='{} {} {}'.format(*origin),
                        SCORING=self.scoring, RADIUS=self.radius)
        inputfile = _TEMPLATE.safe_substitute(replaces)
        inputfilepath = self.tmpfile + '.conf'
        with open(self.tmpfile + '.conf', 'w') as f:
            f.write(inputfile)
        self._paths['conf'] = inputfilepath
        return [self.executable, inputfilepath]

    def parse_output(self, filename):
        """ Get last word of first line (and unique) and parse it into float """
        fitness = 4
        with open(filename) as f:
            for line in f:
                if not line.strip():
                    continue
                fields = line.split()
                if fields[0] == 'Status':
                    fitness = fields.index(self.score_component)
                elif fields[0] == 'Ok':
                    return float(fields[fitness])

    def clean(self):
        for p in self._paths.values():
            os.remove(p)
        self._paths.clear()