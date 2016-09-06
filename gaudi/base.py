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
Contains the core classes we use to build individuals
(potential solutions of the optimization process).

"""

# Python
from collections import OrderedDict
from contextlib import contextmanager
from copy import deepcopy
from zipfile import ZipFile, ZIP_DEFLATED, ZIP_STORED
import logging
import os
import pprint
import sys
# Chimera
import chimera
# External dependencies
import deap.base
import yaml
import numpy as np
# GAUDI
import gaudi.plugin
import gaudi.similarity

pp = pprint.PrettyPrinter(4)

logger = logging.getLogger(__name__)


class Individual(object):

    """
    Base class for `individual` objects that are evaluated by DEAP.

    Each individual is a potential solution. It contains all that is needed
    for an evaluation. With multiprocessing in mind, individuals should be
    self-contained so it can be passed between threads.

    The defined methods are only wrapper calls to the respective methods of each
    gene.

    Parameters
    ----------

    cfg : gaudi.parse.Settings
        The full parsed object from the configuration YAML file.
    cache : dict or dict-like
        A mutable object that can be used to store values across instances.
    dummy : bool
        If True, create an uninitialized Individual, only containing the
        cfg attribute. If false, call `__ready__` and complete initialization.

    Atttributes
    -----------
    __CACHE : dict
        Class attribute that caches gene data across instances
    __CACHE_OBJ : dict
        Class attribute taht caches objectives data across instances

    .. todo::

        :meth:`write` should use Pickle and just save the whole object, but
        Chimera's inmutable objects (Atoms, Residues, etc) get in the way. A
        workaround may be found if we take a look a the session saving code.

    """

    _CACHE = {}
    _CACHE_OBJ = {}

    def __init__(self, cfg=None, cache=None, **kwargs):
        logger.debug("Creating new individual with id %s", id(self))
        self.cfg = cfg
        self.genes = OrderedDict()
        self.fitness = None
        self.expressed = False
        self._molecules = OrderedDict()
        self.__ready__()
        self.__expression_hooks__()
        
    def __ready__(self):
        """
        A *post-init* method used to avoid initialization problems with
        `__deepcopy__`. It's just the second part of a two-stage
        `__init__`.
        """
        if self.cfg is not None:
            gaudi.plugin.load_plugins(self.cfg.genes, container=self.genes,
                                      parent=self,
                                      cxeta=self.cfg.ga.cx_eta,
                                      mteta=self.cfg.ga.cx_eta,
                                      indpb=self.cfg.ga.mut_indpb)

            self.fitness = Fitness(self.cfg.weights)
            mod, fn = self.cfg.similarity.module.rsplit('.', 1)
            self._similarity = getattr(sys.modules[mod], fn)

        for name, gene in self.genes.items():
            gene.__ready__()
            if gene.__class__.__name__ == 'Molecule':
                self._molecules[name] = gene

    def __expression_hooks__(self):
        with expressed(self):
            for name, gene in self.genes.items():
                gene.__expression_hooks__()

    def __deepcopy__(self, memo):
        new = self.__class__(cfg=self.cfg, dummy=True)
        new.genes = deepcopy(self.genes, memo)
        new.fitness = deepcopy(self.fitness, memo)
        new._similarity = self._similarity
        new.expressed = self.expressed
        for name, gene in new.genes.items():
            gene.parent = new
            if gene.__class__.__name__ == 'Molecule':
                new._molecules[name] = gene

        return new

    def evaluate(self, environment):
        """
        Express individual, evaluate it and unexpress it.

        Parameters
        ----------
        environment : Environment
            Objectives that will evaluate the individual
        """
        logger.debug("Evaluating individual #%s", id(self))
        self.express()
        score = environment.evaluate(self)
        self.unexpress()
        return score

    def express(self):
        """
        Express genes in this environment. Very much like 'compiling' the
        individual to a chimera.Molecule.
        """
        for name, gene in self.genes.items():
            logger.debug("Expressing gene %s with allele\n%s",
                         name, pp.pformat(gene.allele))
            gene.express()

        for molecule in self._molecules.values():
            coords = np.array([a.xformCoord() for a in sorted(molecule.compound.mol.atoms,
                                                             key=lambda a: a.serialNumber)])
            molecule._expressed_xformcoords_cache = coords
        self.expressed = True

    def unexpress(self):
        """
        Undo .express()
        """
        for gene in reversed(self.genes.values()):
            logger.debug("Reverting expression of gene %s", gene.name)
            gene.unexpress()
        self.expressed = False

    def mate(self, other):
        """
        Recombine genes of `self` with `other`. It simply calls `mate` on
        each gene instance

        Parameters
        ----------
        other : Individual
            Another individual to mate with.
        """
        for gene in self.genes.values():
            gene.mate(other.genes[gene.name])
        logger.debug("#%s mated #%s", id(self), id(other))
        return self, other

    def mutate(self, indpb):
        """
        Trigger a round of possible mutations across all genes

        Parameters
        ----------
        indpb : float
            Probability of suffering a mutation
        """
        for gene in self.genes.values():
            gene.mutate(indpb)
        logger.debug("#%s mutated", id(self))
        return self,

    def similar(self, other):
        """
        Compare `self` and `other` with a similarity function.

        Returns
        -------
        bool
        """
        return self._similarity(self, other,
                                *self.cfg.similarity.args,
                                **self.cfg.similarity.kwargs)

    def write(self, i):
        """
        Export the individual to a mol2 file

        Parameters
        ----------
        i : int
            Individual identificator in current generation or hall of fame

        .. note :: 

            Maybe someday we can pickle it all :/
            >>> filename = os.path.join(path, '{}_{}.pickle.gz'.format(name,i))
            >>> with gzip.GzipFile(filename, 'wb') as f:
            >>>     cPickle.dump(self, f, 0)
            >>> return filename
        """
        path = self.cfg.output.path
        name = self.cfg.output.name
        COMPRESS = ZIP_DEFLATED if self.cfg.output.compress else ZIP_STORED
        self.express()
        zipfilename = os.path.join(path, '{}_{:03d}.zip'.format(name, i))
        with ZipFile(zipfilename, 'w', COMPRESS) as z:
            output = OrderedDict()
            for gene in self.genes.values():
                logger.debug("Writing %s to file", gene.name)
                filename = gene.write(path, "{}_{:03d}".format(name, i))
                if filename:
                    z.write(filename, os.path.basename(filename))
                    os.remove(filename)
                    output[gene.name] = os.path.basename(filename)
            try:
                output['score'] = list(self.fitness.values)
            except AttributeError:  # fitness not in individual :/
                raise
            z.writestr('{}_{:03d}.gaudi'.format(name, i),
                       yaml.dump(output, default_flow_style=False))
        self.unexpress()
        return zipfilename

    def clear_cache(self):
        self._CACHE.clear()
        self._CACHE_OBJ.clear()
        for gene in self.genes.values():
            gene.clear_cache()


@contextmanager
def expressed(individual):
    individual.express()
    try:
        yield individual
    finally:
        individual.unexpress()


class Environment(object):

    """
    Objective container and helper to evaluate an individual. It must be
    instantiated with a gaudi.parse.Settings object.

    Parameters
    ----------
    cfg : gaudi.parse.Settings
        The parsed configuration YAML file that contains objectives information
    """

    def __init__(self, cfg=None, *args, **kwargs):
        self.zone = chimera.selection.ItemizedSelection()
        self.objectives = OrderedDict()
        self.cfg = cfg
        if self.cfg is not None:
            self.weights = self.cfg.weights
            gaudi.plugin.load_plugins(self.cfg.objectives,
                                      container=self.objectives,
                                      zone=self.zone, environment=self)

    def evaluate(self, individual):
        """
        individual : Individual
        """
        scores = []
        individual.express()
        for name, obj in self.objectives.items():
            score = obj.evaluate(individual)
            scores.append(score)
            logger.debug("%s fitness is %s", name, score)
        individual.unexpress()
        return scores

    def clear_cache(self):
        for obj in self.objectives.values():
            obj.clear_cache()


class Fitness(deap.base.Fitness):

    wvalues = ()

    def __init__(self, weights):
        self.weights = weights
        deap.base.Fitness.__init__(self)

    def __deepcopy__(self, memo):
        new = self.__class__(self.weights)
        new.wvalues = self.wvalues + ()
        return new
