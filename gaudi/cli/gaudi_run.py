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
`gaudi.cli.gaudi_run` is the main hub for launching GAUDI jobs.

It sets up the configuration environment needed by DEAP (responsible for the GA)
and ties it up to the GAUDI custom classes that shape up the invididuals and
objectives. All in a loosely-coupled approach based on Python modules called on-demand.

**Usage**. Simply, type:

.. code-block :: console

    gaudi run /path/to/job.gaudi-input

If the previous does not work, try with the manual mode:

.. code-block :: console

    cd /path/to/gaudi/installation/directory/
    /path/to/chimera/bin/chimera --nogui --script "gaudi_run.py /path/to/job.gaudi-input"


.. todo::

    DEAP default algorithm is nice, but we will be needing some custom features. For
    example, handle KeyboardInterrupt not to lose the population, and so on.

"""

# Python
from __future__ import print_function
from importlib import import_module
import logging
import os
import shutil
import sys
# External dependencies
try:
    import chimera
except ImportError:
    print("Chimera not importable from this environment. Please, install "
          "PyChimera and use it to run this file.")
import numpy
import deap.creator
import deap.tools
import deap.base
import deap.algorithms
import yaml
# from multiprocess import Pool
# GAUDI
import gaudi.algorithms
import gaudi.base
import gaudi.box
import gaudi.genes
import gaudi.objectives
import gaudi.parse
import gaudi.plugin
import gaudi.similarity

if sys.version_info.major == 3:
    basestring = str


def launch(cfg):
    """
    Runs a GAUDI job

    Parameters
    ----------
    cfg : gaudi.parse.Settings
        Parsed YAML dict with attribute-like access
    """
    gaudi.plugin.import_plugins(*cfg.genes)
    gaudi.plugin.import_plugins(*cfg.objectives)
    import_module(cfg.similarity.module.rsplit('.', 1)[0])

    # DEAP setup: Fitness, Individuals, Population
    toolbox = deap.base.Toolbox()
    toolbox.register("call", (lambda fn, *args, **kwargs: fn(*args, **kwargs)))
    toolbox.register("individual", toolbox.call, gaudi.base.MolecularIndividual, cfg)
    toolbox.register("population", deap.tools.initRepeat, list, toolbox.individual)
    population = toolbox.population(n=cfg.ga.population)

    environment = gaudi.base.Environment(cfg)

    toolbox.register("evaluate", lambda ind: environment.evaluate(ind))
    toolbox.register("mate", (lambda ind1, ind2: ind1.mate(ind2)))
    toolbox.register("mutate", (lambda ind, indpb: ind.mutate(indpb)), indpb=cfg.ga.mut_indpb)
    toolbox.register("similarity", (lambda ind1, ind2: ind1.similar(ind2)))
    toolbox.register("select", deap.tools.selNSGA2)
    
    # Multiprocessing pool (Resolve pickle!)
    # pool = Pool(initializer=lambda: setattr(sys, 'stdout', open(str(os.getpid()) + ".out", "w")))
    # toolbox.register("map", pool.map)

    if cfg.output.history:
        history = deap.tools.History()
        # Decorate the variation operators
        toolbox.decorate("mate", history.decorator)
        toolbox.decorate("mutate", history.decorator)
        history.update(population)

    if cfg.output.pareto:
        elite = deap.tools.ParetoFront(toolbox.similarity)
    else:
        elite_size = int(cfg.ga.population * cfg.ga.mu)
        elite = deap.tools.HallOfFame(elite_size, similar=toolbox.similarity)
    
    if cfg.output.verbose:
        stats = deap.tools.Statistics(lambda ind: ind.fitness.values)
        numpy.set_printoptions(precision=cfg.output.precision)
        stats.register("avg", numpy.mean, axis=0)
        stats.register("std", numpy.std, axis=0)
        stats.register("min", numpy.min, axis=0)
        stats.register("max", numpy.max, axis=0)
    else:
        stats = None

    # Begin evolution
    mu = int(cfg.ga.mu * cfg.ga.population)
    lambda_ = int(cfg.ga.lambda_ * cfg.ga.population)
    population, log = gaudi.algorithms.ea_mu_plus_lambda(
        population, toolbox, cfg=cfg, mu=mu, lambda_=lambda_,
        cxpb=cfg.ga.cx_pb, mutpb=cfg.ga.mut_pb, 
        ngen=cfg.ga.generations, halloffame=elite,
        verbose=cfg.output.verbose, stats=stats)

    return population, log, elite


def enable_logging(path=None, name=None, debug=False):
    """
    Register loggers and handlers for both stdout and file
    """

    class CustomFormatter(logging.Formatter):

        CUSTOM_FORMATS = {
            logging.DEBUG: "DEBUG: %(module)s: %(lineno)d: %(message)s",
            logging.INFO: "%(message)s",
            logging.WARNING: "Warning: %(message)s",
            logging.ERROR: "[!] %(message)s",
            logging.CRITICAL: "CRITICAL: %(message)s",
            100: "%(message)s"
        }

        def format(self, record):
            format_orig = self._fmt
            self._fmt = self.CUSTOM_FORMATS.get(record.levelno, format_orig)
            result = logging.Formatter.format(self, record)
            self._fmt = format_orig
            return result

    logger = logging.getLogger('gaudi')
    logger.setLevel(logging.DEBUG)

    # create CONSOLE handler and set level to error
    handler = logging.StreamHandler()
    handler.setLevel(logging.ERROR)
    formatter = CustomFormatter()
    handler.setFormatter(formatter)
    logger.addHandler(handler)

    # create debug file handler and set level to debug
    if path and name:
        handler = logging.FileHandler(os.path.join(path, name + ".gaudi-log"), 'w')
        if debug:
            handler.setLevel(logging.DEBUG)
        else:
            handler.setLevel(15)
        formatter = logging.Formatter(
            "%(asctime)s - %(name)s - %(levelname)s - %(message)s", "%Y.%m.%d %H:%M:%S")
        handler.setFormatter(formatter)
        logger.addHandler(handler)

    return logger


def unbuffer_stdout():
    class Unbuffered(object):

        def __init__(self, stream):
            self.stream = stream

        def write(self, data):
            self.stream.write(data)
            self.stream.flush()

        def __getattr__(self, attr):
            return getattr(self.stream, attr)

    sys.stdout = Unbuffered(sys.stdout)


#@gaudi.box.do_cprofile
def main(cfg, debug=False):
    """
    Starts a GAUDI job

    Parameters
    ----------
    cfg : str or gaudi.parse.Settings
        Path to YAML input file or an already parsed YAML file
        via gaudi.parse.Settings class
    debug : bool, optional, default=False
        Whether to enable verbose logging or not.
    """
    # Parse input file
    if isinstance(cfg, basestring) and os.path.isfile(cfg):
        cfg = gaudi.parse.Settings(cfg)
    
    # Enable logging to stdout and file
    unbuffer_stdout()
    logger = enable_logging(cfg.output.path, cfg.output.name, debug=debug)
    logger.log(100, 'Loaded input %s', cfg._path)

    # Place a copy of input inside output path
    shutil.copy(cfg._path, os.path.join(cfg.output.path, os.path.basename(cfg._path)))

    # Disable Chimera's auto ksdssp
    chimera.triggers.addHandler("Model", gaudi.box.suppress_ksdssp, None)

    # Run simulation
    try:
        logger.log(100, 'Launching job with...')
        logger.log(100, '  Genes: %s', ', '.join([g.name for g in cfg.genes]))
        logger.log(100, '  Objectives: %s', ', '.join([o.name for o in cfg.objectives]))
        pop, log, best = launch(cfg)
    except Exception as e:
        log_path = os.path.join(cfg.output.path, cfg.output.name + ".gaudi-log")
        logger.error('An exception occurred: %s\n    '
                     'Check traceback in logfile %s', e, log_path)
        logger.log(15, "An exception occurred", exc_info=True)
        sys.exit(1)
    gaudi.algorithms.dump_population(best, cfg)

