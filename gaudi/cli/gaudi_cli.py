#!/usr/bin/env python
# -*- coding: utf-8 -*-

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

"""
``gaudi.cli.gaudi_cli`` is the CLI entry point for all GAUDI scripts.

Available commands:

    - run
    - view
    - prepare
    - benchmark
    - rmsd

"""
# Python
import sys
import time
from datetime import timedelta
from importlib import import_module
from textwrap import dedent
# 3rd party
import click
# insilichem
import pychimera
import gaudi


# Helpers
def test_import(command, module):
    try:
        imported = import_module('gaudi.cli.' + module)
    except ImportError as e:
        sys.exit("ERROR: {} (and its dependencies) must be installed"
                 " to use <{}> command.\n{}".format(module, command, e))
    else:
        return imported


def timeit(func, *args, **kwargs):
    def wrapped():
        ts = time.time()
        result = func(*args, **kwargs)
        te = time.time()
        click.echo(
            'Finished after {:0>8}'.format(timedelta(seconds=int(te-ts))))
        return result
    return wrapped


@click.group(context_settings=dict(help_option_names=['-h', '--help']))
@click.version_option(version=gaudi.__version__)
def cli(prog_name='gaudi'):
    """
    GaudiMM: Genetic Algorithms with Unrestricted
    Descriptors for Intuitive Molecular Modeling

    \b
    (C) 2017, InsiliChem
    https://bitbucket.org/insilichem/gaudi
    """
    pychimera.patch_environ()
    pychimera.enable_chimera()
    banner = dedent('''
      .g8"""bgd       db   `7MMF'   `7MF'`7MM"""Yb. `7MMF'                                    
    .dP'     `M      ;MM:    MM       M    MM    `Yb. MM                                      
    dM'       `     ,V^MM.   MM       M    MM     `Mb MM  `7MMpMMMb.pMMMb.  `7MMpMMMb.pMMMb.  
    MM             ,M  `MM   MM       M    MM      MM MM    MM    MM    MM    MM    MM    MM  
    MM.    `7MMF'  AbmmmqMA  MM       M    MM     ,MP MM    MM    MM    MM    MM    MM    MM  
    `Mb.     MM   A'     VML YM.     ,M    MM    ,dP' MM    MM    MM    MM    MM    MM    MM  
      `"bmmmdPY .AMA.   .AMMA.`bmmmmd"'  .JMMmmmdP' .JMML..JMML  JMML  JMML..JMML  JMML  JMML.''')
    click.echo("{}\n{}".format(banner, '-'*len(banner.splitlines()[1])))
    click.echo('GaudiMM: Genetic Algorithms with Unrestricted '
               'Descriptors for Intuitive Molecular Modeling')
    click.echo('{} · v{}\n'.format(gaudi.__copyright__, gaudi.__version__))


@cli.command()
@click.option('--debug', help='Dump debug info to logfile', is_flag=True)
@click.argument('filename', required=True, type=click.Path(exists=True))
def run(filename, debug):
    """
    Launch a GAUDI input file.
    """
    gaudi_run = test_import('run', 'gaudi_run')
    ts = time.time()
    gaudi_run.main(filename, debug)
    te = time.time()
    click.echo('Finished after {:0>8}'.format(timedelta(seconds=int(te-ts))))


@cli.command()
@click.argument('filename', required=False, type=click.Path(exists=True))
def view(filename):
    """
    Analyze the results in a GAUDI output file.
    """
    # gaudinspect = test_import('view', 'gaudinspect')
    # gaudinspect.cli.open_file(filename)
    click.echo("This argument is still unimplemented.")


@cli.command()
@click.argument('filename', required=False, type=click.Path())
def prepare(filename):
    """
    Create or edit a GAUDI input file.
    """
    # gaudinspect = test_import('prepare', 'gaudinspect')
    # gaudinspect.cli.prepare_input(filename)
    click.echo("This argument is still unimplemented.")


@cli.command()
@click.argument('dataset')
@click.argument('template')
def benchmark(dataset, template):
    """
    Performs the same essay over a dataset.
    """
    gaudi_benchmark = test_import('benchmark', 'gaudi_benchmark')
    ts = time.time()
    gaudi_benchmark.main(dataset, template)
    te = time.time()
    click.echo('Finished after {:0>8}'.format(timedelta(seconds=int(te-ts))))


@cli.command()
@click.argument('dataset')
@click.argument('outputfile')
@click.option('--reference', '-r', help='Filename of the reference molecule',
              default='reference.mol2')
@click.option('--results', '-R', help='Subdirectory where results are generated',
              default='results')
@click.option('--great', '-G', help='RMSD threshold in A for a solution to be considered great',
              type=float, default=1.5)
@click.option('--good', '-g', help='RMSD threshold in A for a solution to be considered good',
              type=float, default=2.5)
def rmsd(dataset, outputfile, reference, results, great, good):
    """
    Calculate RMSD of results vs reference.
    """
    gaudi_rmsd = test_import('rmsd', 'gaudi_rmsd')
    ts = time.time()
    gaudi_rmsd.rmsd(dataset, outputfile, reference, results)
    gaudi_rmsd.stats(dataset, great, good)
    te = time.time()
    click.echo('Finished after {:0>8}'.format(timedelta(seconds=(te-ts))))

if "__main__" == __name__:
    cli(prog_name='gaudi')
