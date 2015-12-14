#!/usr/bin/env python
# -*- coding: utf-8 -*-

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
# 3rd party
import click
# GAUDI
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
    GAUDI: Genetic Algorithms for Universal Design Inference

    \b
    By Jaime Rodríguez-Guerra and Jean-Didier Maréchal.
    https://bitbucket.org/jrgp/gaudi
    """
    click.echo('\n')


@cli.command()
@click.argument('filename', required=True, type=click.Path(exists=True))
def run(filename):
    """
    Launch a GAUDI input file.
    """
    gaudi_run = test_import('run', 'gaudi_run')
    ts = time.time()
    gaudi_run.main(filename)
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
@click.option('--gaudi', '-g', help='Location of gaudi binary',
              default='gaudi')
def benchmark(dataset, template, gaudi):
    """
    Performs the same essay over a dataset.
    """
    gaudi_benchmark = test_import('benchmark', 'gaudi_benchmark')
    ts = time.time()
    gaudi_benchmark.main(dataset, template, gaudi)
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
