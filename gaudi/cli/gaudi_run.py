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
`gaudi.cli.gaudi_run` is the main hub for launching GAUDI essays.

It sets up the configuration environment needed by DEAP (responsible for the GA)
and ties it up to the GAUDI custom classes that shape up the invididuals and
objectives. All in a loosely-coupled approach based on Python modules called on-demand.

**Usage**. Simply, type:

.. code-block :: console

    gaudi run /path/to/essay.gaudi-input

If the previous does not work, try with the manual mode:

.. code-block :: console

    cd /path/to/gaudi/installation/directory/
    /path/to/chimera/bin/chimera --nogui --script "gaudi_run.py /path/to/essay.gaudi-input"


.. todo::

    DEAP default algorithm is nice, but we will be needing some custom features. For
    example, handle KeyboardInterrupt not to lose the population, and so on.

"""

# Python
from __future__ import print_function
from time import strftime
from importlib import import_module
import logging
import numpy
import os
import sys
# External dependencies
try:
    import chimera
except ImportError:
    print("You must install UCSF Chimera and run GAUDI with its own Python interpreter.\n"
          "Check the install guide for more details.")
import deap.creator
import deap.tools
import deap.base
import deap.algorithms
import yaml
# GAUDI
import gaudi.algorithms
import gaudi.base
import gaudi.box
import gaudi.genes
import gaudi.objectives
import gaudi.parse
import gaudi.plugin
import gaudi.similarity
import gaudi.version


def launch(cfg):
    """
    Runs a GAUDI essay

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
    toolbox.register("individual", toolbox.call, gaudi.base.Individual, cfg)
    toolbox.register("population", deap.tools.initRepeat, list, toolbox.individual)
    population = toolbox.population(n=cfg.ga.population)

    environment = gaudi.base.Environment(cfg)

    toolbox.register("evaluate", lambda ind: environment.evaluate(ind))
    toolbox.register("mate", (lambda ind1, ind2: ind1.mate(ind2)))
    toolbox.register("mutate", (lambda ind, indpb: ind.mutate(indpb)), indpb=cfg.ga.mut_indpb)
    toolbox.register("similarity", (lambda ind1, ind2: ind1.similar(ind2)))
    toolbox.register("select", deap.tools.selNSGA2)

    if cfg.output.history:
        history = deap.tools.History()
        # Decorate the variation operators
        toolbox.decorate("mate", history.decorator)
        toolbox.decorate("mutate", history.decorator)
        history.update(population)

    if cfg.output.pareto:
        best_individuals = deap.tools.ParetoFront(toolbox.similarity)
    else:
        # hof_size_percent = int(0.1 * cfg.ga.population)
        # hof_size = hof_size_percent if hof_size_percent > 2 else 2
        hof_size = int(cfg.ga.population * cfg.ga.mu)
        best_individuals = deap.tools.HallOfFame(hof_size, similar=toolbox.similarity)
    stats = deap.tools.Statistics(lambda ind: ind.fitness.values)
    numpy.set_printoptions(precision=cfg.output.precision)
    stats.register("avg", numpy.mean, axis=0)
    stats.register("min", numpy.min, axis=0)
    stats.register("max", numpy.max, axis=0)

    # Begin evolution
    population, log = deap.algorithms.eaMuPlusLambda(
        population, toolbox, mu=int(cfg.ga.mu * cfg.ga.population),
        lambda_=int(cfg.ga.lambda_ * cfg.ga.population), cxpb=cfg.ga.cx_pb, mutpb=cfg.ga.mut_pb,
        ngen=cfg.ga.generations, stats=stats, halloffame=best_individuals)

    return population, log, best_individuals



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
        handler = logging.FileHandler(os.path.join(path, name + ".gaudi-log"), "w")
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
def main(filename, debug=False):
    # Parse input file
    cfg = gaudi.parse.Settings(filename)
    # Enable logging to stdout and file
    unbuffer_stdout()
    logger = enable_logging(cfg.output.path, cfg.output.name, debug=debug)
    logger.log(100, 'Loaded input %s', filename)

    # Disable auto ksdssp
    chimera.triggers.addHandler("Model", gaudi.box.suppress_ksdssp, None)

    # Run simulation
    try:
        logger.log(100, 'Launching essay ...')
        pop, log, best = launch(cfg)
    except Exception as e:
        log_path = os.path.join(cfg.output.path, cfg.output.name + ".gaudi-log")
        logger.error('An exception occurred: %s\n    '
                     'Check traceback in logfile %s', e, log_path)
        logger.log(15, "An exception occurred", exc_info=True)
        sys.exit(1)

    # Write results
    logger.log(100, 'Writing %s results to disk', len(pop))
    results = {'GAUDI.objectives': [
        '{} ({})'.format(obj.name, obj.module) for obj in cfg.objectives]}
    results['GAUDI.results'] = {}
    for i, ind in enumerate(best):
        filename = ind.write(i)
        results['GAUDI.results'][os.path.basename(filename)] = map(float, ind.fitness.values)

    outputpath = os.path.join(cfg.output.path, '{}.gaudi-output'.format(cfg.output.name))
    with open(outputpath, 'w+') as out:
        out.write('# Generated by GAUDI on {}\n\n'.format(strftime("%Y-%m-%d %H:%M:%S")))
        out.write(yaml.safe_dump(results, default_flow_style=False))

if __name__ == "__main__":
    try:
        main(sys.argv[1])
    except IndexError:
        sys.exit("ERROR: Input file not provided. \n")
